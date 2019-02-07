# -*- coding: utf-8 -*-

"""Assemble a BEL graph as an `ideogram <https://github.com/eweitz/ideogram>`_ chart in HTML.."""

import os
from typing import Any, Mapping, Set

import pandas as pd
from jinja2 import Template

import bio2bel_hgnc
from bio2bel_hgnc.models import HumanGene
from pybel import BELGraph
from pybel.constants import GENE
from pybel.dsl import Gene
from pybel.struct import get_nodes_by_function
from pybel.struct.mutation import collapse_all_variants, enrich_protein_and_rna_origins

__all__ = [
    'render',
]

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)))
COLUMNS = ['start_position_on_the_genomic_accession',
           'end_position_on_the_genomic_accession', ]

# Create a template object from the template file, load once
TEMPLATE_PATH = os.path.join(HERE, 'index.html')
with open(TEMPLATE_PATH, 'rt') as f:
    template_str = f.read()
    template = Template(template_str)


def render(graph: BELGraph) -> str:
    """Render the graph as an HTML string."""
    annotations = list(prerender(graph).values())
    return template.render(
        annotations=annotations,
        title=graph.name or 'BEL Graph Information Density',
    )


def prerender(graph: BELGraph) -> Mapping[str, Mapping[str, Any]]:
    """Generate the annotations JSON for Ideogram."""
    graph: BELGraph = graph.copy()
    enrich_protein_and_rna_origins(graph)
    collapse_all_variants(graph)
    genes: Set[Gene] = get_nodes_by_function(graph, GENE)
    hgnc_symbols = {
        gene.name
        for gene in genes
        if gene.namespace.lower() == 'hgnc'
    }

    result = {}

    hgnc_manager = bio2bel_hgnc.Manager()
    human_genes = (
        hgnc_manager.session
            .query(HumanGene.symbol, HumanGene.location)
            .filter(HumanGene.symbol.in_(hgnc_symbols))
            .all()
    )
    for human_gene in human_genes:
        result[human_gene.symbol] = {
            'name': human_gene.symbol,
            'chr': (
                human_gene.location.split('q')[0]
                if 'q' in human_gene.location else
                human_gene.location.split('p')[0]
            ),
        }

    df = get_df()

    for _, (gene_id, symbol, start, stop) in df[df['Symbol'].isin(hgnc_symbols)].iterrows():
        result[symbol]['start'] = start
        result[symbol]['stop'] = stop

    return result


def get_df() -> pd.DataFrame:
    return pd.read_csv(os.path.join(HERE, 'refseq_human.csv'))
