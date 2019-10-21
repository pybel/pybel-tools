# -*- coding: utf-8 -*-

"""Assemble a BEL graph as an `ideogram <https://github.com/eweitz/ideogram>`_ chart in HTML.."""

import random
from typing import Any, Mapping, Optional, TextIO

from IPython.display import Javascript

import bio2bel_hgnc
from bio2bel_entrez.parser import get_human_refseq_slim_df
from bio2bel_hgnc.models import HumanGene
from pybel import BELGraph
from pybel.dsl import CentralDogma
from ..jinja_utils import build_template_renderer

__all__ = [
    'to_html',
    'to_html_file',
    'to_html_path',
    'to_jupyter',
]

COLUMNS = [
    'start_position_on_the_genomic_accession',
    'end_position_on_the_genomic_accession',
]


render_template = build_template_renderer(__file__)


def to_jupyter(graph: BELGraph, chart: Optional[str] = None) -> Javascript:
    """Render the graph as JavaScript in a Jupyter Notebook."""
    context = _get_context(graph, chart=chart)
    javascript_str = render_template('render_with_javascript.js', **context)
    return Javascript(javascript_str)


def to_html(graph: BELGraph, chart: Optional[str] = None) -> str:
    """Render the graph as an HTML string.

    Common usage may involve writing to a file like:

    >>> from pybel.examples import sialic_acid_graph
    >>> with open('ideogram_output.html', 'w') as file:
    ...     print(to_html(sialic_acid_graph), file=file)
    """
    context = _get_context(graph, chart=chart)
    return render_template('index.html', **context)


def to_html_file(graph: BELGraph, file: Optional[TextIO] = None, chart: Optional[str] = None) -> None:
    """Write the graph as an HTML file."""
    print(to_html(graph=graph, chart=chart), file=file)


def to_html_path(graph: BELGraph, path: str, chart: Optional[str] = None) -> None:
    """Write the graph as an HTML file."""
    with open(path, 'w') as file:
        to_html_file(graph=graph, file=file, chart=chart)


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


def prerender(graph: BELGraph, hgnc_manager=None) -> Mapping[str, Mapping[str, Any]]:
    """Generate the annotations JSON for Ideogram."""
    if hgnc_manager is None:
        hgnc_manager = bio2bel_hgnc.Manager()

    hgnc_symbols = {
        node.name
        for node in graph
        if isinstance(node, CentralDogma) and node.namespace.lower() == 'hgnc'
    }

    refseq_df = get_human_refseq_slim_df()

    result = {
        hgnc_symbol: dict(name=hgnc_symbol, start=start, stop=stop)
        for _, hgnc_symbol, start, stop in refseq_df[refseq_df['Symbol'].isin(hgnc_symbols)].values
    }

    human_genes = (
        hgnc_manager.session
        .query(HumanGene.symbol, HumanGene.location)
        .filter(HumanGene.symbol.in_(hgnc_symbols))
        .all()
    )
    for human_gene in human_genes:
        if human_gene.symbol not in result:
            continue  # something doesn't have a mapping in HGNC
        result[human_gene.symbol]['chr'] = (
            human_gene.location.split('q')[0]
            if 'q' in human_gene.location else
            human_gene.location.split('p')[0]
        )

    return result
