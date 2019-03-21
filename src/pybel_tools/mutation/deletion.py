# -*- coding: utf-8 -*-

"""Deletion functions to supplement :mod:`pybel.struct.mutation.deletion`."""

from typing import Iterable, Tuple

from pybel import BELGraph
from pybel.dsl import BaseEntity
from pybel.struct.pipeline import in_place_transformation
from ..summary.edge_summary import pair_is_consistent

__all__ = [
    'remove_inconsistent_edges',
]


@in_place_transformation
def remove_inconsistent_edges(graph: BELGraph) -> None:
    """Remove all edges between node pairs with inconsistent edges.

    This is the all-or-nothing approach. It would be better to do more careful investigation of the evidences during
    curation.
    """
    for u, v in get_inconsistent_edges(graph):
        edges = [(u, v, k) for k in graph[u][v]]
        graph.remove_edges_from(edges)


def get_inconsistent_edges(graph: BELGraph) -> Iterable[Tuple[BaseEntity]]:
    """Iterate over pairs of nodes with inconsistent edges."""
    for u, v in graph.edges():
        if not pair_is_consistent(graph, u, v):
            yield u, v
