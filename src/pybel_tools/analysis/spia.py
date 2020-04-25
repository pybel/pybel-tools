# -*- coding: utf-8 -*-

"""An exporter for signaling pathway impact analysis (SPIA) described by [Tarca2009]_.

.. [Tarca2009] Tarca, A. L., *et al* (2009). `A novel signaling pathway impact analysis
               <https://doi.org/10.1093/bioinformatics/btn577>`_. Bioinformatics, 25(1), 75–82.

To run this module on an arbitrary BEL graph, use the command ``python -m pybel_tools.analysis.spia``.

.. seealso:: https://bioconductor.org/packages/release/bioc/html/SPIA.html
"""

import sys

import click

from pybel import BELGraph
from pybel.cli import graph_pickle_argument
from pybel.io.spia import spia_matrices_to_excel, spia_matrices_to_tsvs, to_spia_dfs as bel_to_spia_matrices

__all__ = [
    'bel_to_spia_matrices',
    'spia_matrices_to_excel',
    'spia_matrices_to_tsvs',
    'spia',
]


@click.command()
@graph_pickle_argument
@click.option('--xlsx', type=click.Path(file_okay=True, dir_okay=False))
@click.option('--tsvs', type=click.Path(file_okay=False, dir_okay=True))
def spia(graph: BELGraph, xlsx: str, tsvs: str):
    """Export the graph to a SPIA Excel sheet."""
    if not xlsx and not tsvs:
        click.secho('Specify at least one option --xlsx or --tsvs', fg='red')
        sys.exit(1)

    spia_matrices = bel_to_spia_matrices(graph)

    if xlsx:
        spia_matrices_to_excel(spia_matrices, xlsx)

    if tsvs:
        spia_matrices_to_tsvs(spia_matrices, tsvs)


if __name__ == '__main__':
    spia()
