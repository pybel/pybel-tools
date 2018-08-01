# -*- coding: utf-8 -*-

"""This module contains convenient functions for removing nodes/edges that are returned from selection functions"""

from pybel.struct.pipeline import in_place_transformation

from ..summary.edge_summary import pair_is_consistent

__all__ = [
    'remove_inconsistent_edges',
]


def get_inconsistent_edges(graph):
    """Returns an iterator over inconsistent edges

    :param pybel.BELGraph graph: A BEL graph
    :return: An iterator over (source, target) node pairs corresponding to edges with many inconsistent relations
    :rtype: iter
    """
    for u, v in graph.edges():
        if not pair_is_consistent(graph, u, v):
            yield u, v


@in_place_transformation
def remove_inconsistent_edges(graph):
    """Remove all edges between node paris with consistent edges.

    This is the all-or-nothing approach. It would be better to do more careful investigation of the evidences during
    curation.

    :param pybel.BELGraph graph: A BEL graph
    """
    for u, v in get_inconsistent_edges(graph):
        edges = [(u, v, k) for k in graph[u][v]]
        graph.remove_edges_from(edges)
