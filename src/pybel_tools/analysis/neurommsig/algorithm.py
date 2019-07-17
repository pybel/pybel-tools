# -*- coding: utf-8 -*-

"""An implementation of the NeuroMMSig mechanism enrichment algorithm [DomingoFernandez2017]_.

.. [DomingoFernandez2017] ﻿Domingo-Fernández, D., *et al* (2017). `Multimodal mechanistic signatures for
    neurodegenerative diseases (NeuroMMSig): A web server for mechanism enrichment
    <https://doi.org/10.1093/bioinformatics/btx399>`_. Bioinformatics, 33(22), 3679–3681.
"""

import itertools as itt
import logging
from collections import Counter
from typing import List, Mapping, Optional

from tqdm import tqdm

from pybel import BELGraph, Pipeline
from pybel.constants import GENE
from pybel.dsl import BaseEntity, Gene
from pybel.struct import (
    collapse_all_variants, collapse_to_genes, enrich_protein_and_rna_origins, get_nodes_by_function,
    get_subgraphs_by_annotation,
)
from ...utils import calculate_betweenness_centality

__all__ = [
    'get_neurommsig_scores',
    'get_neurommsig_score',
    'neurommsig_graph_preprocessor',
]

logger = logging.getLogger(__name__)

neurommsig_graph_preprocessor = Pipeline.from_functions([
    enrich_protein_and_rna_origins,
    collapse_to_genes,
    collapse_all_variants,
])


def get_neurommsig_scores(
        graph: BELGraph,
        genes: List[Gene],
        annotation: str = 'Subgraph',
        ora_weight: Optional[float] = None,
        hub_weight: Optional[float] = None,
        top_percent: Optional[float] = None,
        topology_weight: Optional[float] = None,
        preprocess: bool = False,
        use_tqdm: bool = False,
        tqdm_kwargs: Optional[Mapping] = None,
) -> Optional[Mapping[str, float]]:
    """Preprocess the graph, stratify by the given annotation, then run the NeuroMMSig algorithm on each.

    :param graph: A BEL graph
    :param genes: A list of gene nodes
    :param annotation: The annotation to use to stratify the graph to subgraphs
    :param ora_weight: The relative weight of the over-enrichment analysis score from
     :py:func:`neurommsig_gene_ora`. Defaults to 1.0.
    :param hub_weight: The relative weight of the hub analysis score from :py:func:`neurommsig_hubs`.
     Defaults to 1.0.
    :param top_percent: The percentage of top genes to use as hubs. Defaults to 5% (0.05).
    :param topology_weight: The relative weight of the topolgical analysis core from
     :py:func:`neurommsig_topology`. Defaults to 1.0.
    :param preprocess: If true, preprocess the graph.
    :return: A dictionary from {annotation value: NeuroMMSig composite score}

    Pre-processing steps:

    1. Infer the central dogma with :func:``
    2. Collapse all proteins, RNAs and miRNAs to genes with :func:``
    3. Collapse variants to genes with :func:``
    """
    if preprocess:
        graph = neurommsig_graph_preprocessor.run(graph)

    if all(isinstance(gene, str) for gene in genes):
        genes = [Gene('HGNC', gene) for gene in genes]

    if all(gene not in graph for gene in genes):
        logger.warning('no genes mapping to graph')
        return

    subgraphs = get_subgraphs_by_annotation(graph, annotation=annotation)

    return get_neurommsig_scores_prestratified(
        subgraphs=subgraphs,
        genes=genes,
        ora_weight=ora_weight,
        hub_weight=hub_weight,
        top_percent=top_percent,
        topology_weight=topology_weight,
        use_tqdm=use_tqdm,
        tqdm_kwargs=tqdm_kwargs,
    )


def get_neurommsig_scores_prestratified(
        subgraphs: Mapping[str, BELGraph],
        genes: List[Gene],
        ora_weight: Optional[float] = None,
        hub_weight: Optional[float] = None,
        top_percent: Optional[float] = None,
        topology_weight: Optional[float] = None,
        use_tqdm: bool = False,
        tqdm_kwargs: Optional[Mapping] = None,
) -> Mapping[str, float]:
    """Run NeuroMMSig on a graph strata.

    :param subgraphs: A pre-stratified set of graphs
    :param genes: A list of gene nodes
    :param ora_weight: The relative weight of the over-enrichment analysis score from
     :py:func:`neurommsig_gene_ora`. Defaults to 1.0.
    :param hub_weight: The relative weight of the hub analysis score from :py:func:`neurommsig_hubs`.
     Defaults to 1.0.
    :param top_percent: The percentage of top genes to use as hubs. Defaults to 5% (0.05).
    :param topology_weight: The relative weight of the topolgical analysis core from
     :py:func:`neurommsig_topology`. Defaults to 1.0.
    :param use_tqdm: If true, show a progress bar
    :return: A dictionary from {annotation value: NeuroMMSig composite score}

    Pre-processing steps:

    1. Infer the central dogma with :func:``
    2. Collapse all proteins, RNAs and miRNAs to genes with :func:``
    3. Collapse variants to genes with :func:``
    """
    it = subgraphs.items()
    if use_tqdm:
        it = tqdm(it, **(tqdm_kwargs or {}))
    return {
        name: get_neurommsig_score(
            graph=subgraph,
            genes=genes,
            ora_weight=ora_weight,
            hub_weight=hub_weight,
            top_percent=top_percent,
            topology_weight=topology_weight,
        )
        for name, subgraph in it
    }


def get_neurommsig_score(
        graph: BELGraph,
        genes: List[Gene],
        ora_weight: Optional[float] = None,
        hub_weight: Optional[float] = None,
        top_percent: Optional[float] = None,
        topology_weight: Optional[float] = None,
) -> float:
    """Calculate the composite NeuroMMSig Score for a given list of genes.

    :param graph: A BEL graph
    :param genes: A list of gene nodes
    :param ora_weight: The relative weight of the over-enrichment analysis score from
     :py:func:`neurommsig_gene_ora`. Defaults to 1.0.
    :param hub_weight: The relative weight of the hub analysis score from :py:func:`neurommsig_hubs`.
     Defaults to 1.0.
    :param top_percent: The percentage of top genes to use as hubs. Defaults to 5% (0.05).
    :param topology_weight: The relative weight of the topolgical analysis core from
     :py:func:`neurommsig_topology`. Defaults to 1.0.
    :return: The NeuroMMSig composite score
    """
    ora_weight = ora_weight or 1.0
    hub_weight = hub_weight or 1.0
    topology_weight = topology_weight or 1.0
    total_weight = ora_weight + hub_weight + topology_weight

    genes = list(genes)

    ora_score = neurommsig_gene_ora(graph, genes)
    hub_score = neurommsig_hubs(graph, genes, top_percent=top_percent)
    topology_score = neurommsig_topology(graph, genes)

    weighted_sum = (
        ora_weight * ora_score +
        hub_weight * hub_score +
        topology_weight * topology_score
    )

    return weighted_sum / total_weight


def neurommsig_gene_ora(graph: BELGraph, genes: List[Gene]) -> float:
    """Calculate the percentage of target genes mappable to the graph.

    Assume: graph central dogma inferred, collapsed to genes, collapsed variants
    """
    graph_genes = set(get_nodes_by_function(graph, GENE))
    return len(graph_genes.intersection(genes)) / len(graph_genes)


def neurommsig_hubs(graph: BELGraph, genes: List[Gene], top_percent: Optional[float] = None) -> float:
    """Calculate the percentage of target genes mappable to the graph.

    Assume: graph central dogma inferred, collapsed to genes, collapsed variants, graph has more than 20 nodes

    :param graph: A BEL graph
    :param genes: A list of nodes
    :param top_percent: The percentage of top genes to use as hubs. Defaults to 5% (0.05).
    """
    top_percent = top_percent or 0.05

    if graph.number_of_nodes() < 20:
        logger.debug('Graph has less than 20 nodes')
        return 0.0

    graph_genes = set(get_nodes_by_function(graph, GENE))

    bc = Counter({
        node: betweenness_centrality
        for node, betweenness_centrality in calculate_betweenness_centality(graph).items()
        if node in graph_genes
    })

    # TODO consider continuous analog with weighting by percentile
    number_central_nodes = int(len(graph_genes) * top_percent)

    if number_central_nodes < 1:
        number_central_nodes = 1

    number_mappable_central_nodes = sum(
        node in genes
        for node in bc.most_common(number_central_nodes)
    )

    return number_mappable_central_nodes / number_central_nodes


def neurommsig_topology(graph: BELGraph, nodes: List[BaseEntity]) -> float:
    r"""Calculate the node neighbor score for a given list of nodes without considering self-loops.

    .. math::

         \frac{\sum_i^n N_G[i]}{n*(n-1)}
    """
    nodes = list(nodes)
    number_nodes = len(nodes)

    if number_nodes <= 1:
        # log.debug('')
        return 0.0

    unnormalized_sum = sum(
        u in graph[v]
        for u, v in itt.product(nodes, repeat=2)
        if v in graph and u != v
    )

    return unnormalized_sum / (number_nodes * (number_nodes - 1.0))
