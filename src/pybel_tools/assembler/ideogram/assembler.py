# -*- coding: utf-8 -*-

"""Assemble a BEL graph as an `ideogram <https://github.com/eweitz/ideogram>`_ chart in HTML.."""

import os
import random
from typing import Any, Mapping, Optional, Set

import pandas as pd
from IPython.display import Javascript
from jinja2 import Template

from pybel import BELGraph
from pybel.constants import GENE
from pybel.dsl import Gene
from pybel.struct import get_nodes_by_function
from pybel.struct.mutation import collapse_all_variants, enrich_protein_and_rna_origins

__all__ = [
    'to_html',
    'to_jupyter',
]

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)))
COLUMNS = ['start_position_on_the_genomic_accession',
           'end_position_on_the_genomic_accession', ]


def to_jupyter(graph: BELGraph, chart: Optional[str] = None) -> Javascript:
    """Render the graph as JavaScript in a Jupyter Notebook."""
    with open(os.path.join(HERE, 'render_with_javascript.js'), 'rt') as f:
        js_template = Template(f.read())

    return Javascript(js_template.render(**_get_context(graph, chart=chart)))


def to_html(graph: BELGraph, chart: Optional[str] = None) -> str:
    """Render the graph as an HTML string.

    Common usage may involve writing to a file like:

    >>> from pybel.examples import sialic_acid_graph
    >>> with open('ideogram_output.html', 'w') as file:
    ...     print(to_html(sialic_acid_graph), file=file)
    """
    with open(os.path.join(HERE, 'index.html'), 'rt') as f:
        html_template = Template(f.read())

    return html_template.render(**_get_context(graph, chart=chart))


def _get_context(graph: BELGraph, chart: Optional[str] = None) -> Mapping[str, Any]:
    annotations = list(prerender(graph).values())
    return dict(
        annotations=annotations,
        title=graph.name or 'BEL Graph Information Density',
        chart=chart or _generate_id(),
    )


def _generate_id() -> str:
    """Generate a random string of letters."""
    return ''.join(random.sample('abcdefghjkmopqrstuvqxyz', 16))


def prerender(graph: BELGraph) -> Mapping[str, Mapping[str, Any]]:
    """Generate the annotations JSON for Ideogram."""
    import bio2bel_hgnc
    from bio2bel_hgnc.models import HumanGene

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
    """Get RefSeq information as a dataframe."""
    return pd.read_csv(os.path.join(HERE, 'refseq_human.csv'))
