# -*- coding: utf-8 -*-

"""This module has tools for downloading and structuring gene orthology data from HGNC, RGD, and MGI"""

from __future__ import print_function

import pandas as pd
import requests

from pybel.constants import CITATION, CITATION_NAME, CITATION_REFERENCE, CITATION_TYPE, EVIDENCE, ANNOTATIONS
from pybel.constants import GENE, ORTHOLOGOUS, RELATION
from pybel.struct.filters import filter_edges
from . import pipeline
from .constants import PUBMED
from .filters.edge_filters import build_relation_filter
from .utils import safe_add_edge

HGNC = 'HGNC'
MGI = 'MGI'
RGD = 'RGD'

#: Annotations from HGNC
#: Columns: HGNC ID, HGNC Symbol, MGI Curated, MGI Dump, RGD Dump
FULL_RESOURCE = 'http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_mgd_id&col=md_mgd_id&col=md_rgd_id&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&submit=submit'

#: Annotations from HGNC
#: Columns: HGNC Symbol, MGI Symbols
MGI_ONLY = 'http://www.genenames.org/cgi-bin/download?col=gd_app_sym&col=md_mgd_id&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&submit=submit'

#: Annotations from HGNC
#: Columns: HGNC Symbol, RGD Symbols
RGD_ONLY = 'http://www.genenames.org/cgi-bin/download?col=gd_app_sym&col=md_rgd_id&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&submit=submit'

#: Annotations from the Jackson Lab, that maintains the Mouse Genome Informatics Database
MGI_ANNOTATIONS = 'http://www.informatics.jax.org/downloads/reports/MGI_MRK_Coord.rpt'

#: Columns: Human Marker Symbol, Human Entrez Gene ID, HomoloGene ID, Mouse Marker Symbol, MGI Marker Accession ID, High-level Mammalian Phenotype ID (space-delimited)
#: See: http://www.informatics.jax.org/downloads/reports/index.html#pheno
MGI_ORTHOLOGY = 'http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt'

#: Columns: RAT_GENE_SYMBOL	RAT_GENE_RGD_ID	RAT_GENE_NCBI_GENE_ID	HUMAN_ORTHOLOG_SYMBOL	HUMAN_ORTHOLOG_RGD	HUMAN_ORTHOLOG_NCBI_GENE_ID	HUMAN_ORTHOLOG_SOURCE	MOUSE_ORTHOLOG_SYMBOL	MOUSE_ORTHOLOG_RGD	MOUSE_ORTHOLOG_NCBI_GENE_ID	MOUSE_ORTHOLOG_MGI	MOUSE_ORTHOLOG_SOURCE	HUMAN_ORTHOLOG_HGNC_ID
#: First 52 rows are comments with # at beginning and line 53 is the header
RGD_ORTHOLOGY = 'ftp://ftp.rgd.mcw.edu/pub/data_release/RGD_ORTHOLOGS.txt'


def download_orthologies_from_hgnc(path):
    """Downloads the full dump to the given path

    :param path: output path
    """
    res = requests.get(FULL_RESOURCE)

    with open(path, 'w') as file:
        for line in res.iter_lines(decode_unicode=True):
            print(line, file=file)


def structure_orthologies_from_hgnc(lines=None):
    """Structures the orthology data to two lists of pairs of (HGNC, MGI) and (HGNC, RGD) identifiers

    :param lines: The iterable over the downloaded orthologies from HGNC. If None, downloads from HGNC
    :return:
    """
    if lines is None:
        lines = requests.get(FULL_RESOURCE).iter_lines(decode_unicode=True)

    mgi_orthologies = []
    rgd_orthologies = []

    for line in lines:
        hgnc_id, hgnc_symbol, _, mgis, rgds = line.strip().split('\t')

        for mgi in mgis.split(','):
            mgi = mgi.strip()
            mgi = mgi.replace('MGI:', '')
            mgi_orthologies.append((hgnc_symbol, mgi))

        for rgd in rgds.split(','):
            rgd = rgd.strip()
            rgd = rgd.replace('RGD:', '')
            rgd_orthologies.append((hgnc_symbol, rgd))

    return mgi_orthologies, rgd_orthologies


def structure_orthologies_from_rgd(path=None):
    df = pd.read_csv(RGD_ORTHOLOGY if path is None else path, skiprows=52, sep='\t')

    mgi_orthologies = []
    rgd_orthologies = []

    for _, hgnc, rat, mouse in df[['HUMAN_ORTHOLOG_SYMBOL', 'RAT_GENE_SYMBOL', 'MOUSE_ORTHOLOG_SYMBOL']].itertuples():
        mgi_orthologies.append((hgnc, mouse))
        rgd_orthologies.append((hgnc, rat))

    return mgi_orthologies, rgd_orthologies


@pipeline.in_place_mutator
def add_orthology_statements(graph, orthologies, namespace):
    """Adds orthology statements for all orthologous nodes to HGNC nodes

    :param graph:
    :type graph: pybel.BELGraph
    :param orthologies: An iterable over pairs of (HGNC, ORTHOLOG) identifiers
    :type orthologies: list
    """
    for hgnc, ortholog in orthologies:
        hgnc_node = GENE, HGNC, hgnc
        ortholog_node = GENE, namespace, ortholog

        if ortholog_node not in graph:
            continue

        if hgnc_node not in graph:
            graph.add_simple_node(*hgnc_node)

        graph.add_edge(hgnc_node, ortholog_node, attr_dict={
            RELATION: ORTHOLOGOUS,
            CITATION: {
                CITATION_TYPE: PUBMED,
                CITATION_REFERENCE: '25355511',
                CITATION_NAME: 'Rat Genome Database'
            },
            EVIDENCE: 'Asserted from: {}'.format(RGD_ORTHOLOGY),
            ANNOTATIONS: {}
        })


@pipeline.in_place_mutator
def add_mgi_orthology_statements(graph, mgi_orthologies):
    add_orthology_statements(graph, mgi_orthologies, MGI)


@pipeline.in_place_mutator
def add_rgd_orthology_statements(graph, rgd_orthologies):
    add_orthology_statements(graph, rgd_orthologies, RGD)


@pipeline.in_place_mutator
def integrate_orthologies_from_hgnc(graph, lines=None):
    """Adds orthology statements to graph using HGNC symbols, MGI IDs, and RGD IDs.

    For MGI symbols and RGD symbols, use :func:`integrate_orthologies_from_rgd`

    :param graph: A BEL Graph
    :type graph: pybel.BELGraph
    :param lines:
    """
    mgio, rgdo = structure_orthologies_from_hgnc(lines=lines)
    add_mgi_orthology_statements(graph, mgio)
    add_rgd_orthology_statements(graph, rgdo)


@pipeline.in_place_mutator
def integrate_orthologies_from_rgd(graph, path=None):
    """Adds orthology statements to graph using HGNC symbols, MGI symbols, and RGD symbols.

    For MGI IDs and RGD IDs, use :func:`integrate_orthologies_from_hgnc`

    :param graph: A BEL Graph
    :type graph: pybel.BELGraph
    :param path: optional path to local RGD_ORTHOLOGS.txt.
                 Defaults to downloading directly from RGD FTP server with pandas
    """
    mgio, rgdo = structure_orthologies_from_rgd(path=path)
    add_mgi_orthology_statements(graph, mgio)
    add_rgd_orthology_statements(graph, rgdo)


ortholog_filter = build_relation_filter(ORTHOLOGOUS)


@pipeline.in_place_mutator
def collapse_orthologies(graph):
    """Collapses all orthology relations.

    Assumes: orthologies are annotated for edge (u,v) where u is the higher priority node

    :param graph: A BEL Graph
    :type graph: pybel.BELGraph

    .. warning:: This won't work for two way orthology annotations, so it's best to use :func:`integrate_orthologies_from_rgd` first
    """

    orthologs = []

    # TODO rewrite with collapsing functions
    for hgnc, ortholog, _, _ in list(filter_edges(graph, ortholog_filter)):

        for u, _, k, d in graph.in_edges_iter(ortholog, data=True, keys=True):
            if d[RELATION] == ORTHOLOGOUS:
                continue
            safe_add_edge(graph, u, hgnc, k, d)

        for _, v, k, d in graph.out_edges_iter(ortholog, data=True, keys=True):
            safe_add_edge(graph, hgnc, v, k, d)

        orthologs.append(ortholog)

    graph.remove_nodes_from(orthologs)
