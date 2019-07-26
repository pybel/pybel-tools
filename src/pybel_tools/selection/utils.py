# -*- coding: utf-8 -*-

"""Utility gunctions for :mod:`pybel_tools.selection`."""

from typing import Iterable, Optional

from pybel import BELGraph, BaseEntity

__all__ = [
    'get_leaves_by_type',
]


def get_leaves_by_type(
        graph: BELGraph,
        func: Optional[str] = None,
        prune_threshold: int = 1,
) -> Iterable[BaseEntity]:
    """Iterate over all nodes in graph (in-place) with only a connection to one node.

    Useful for gene and RNA. Allows for optional filter by function type.

    :param pybel.BELGraph graph: A BEL graph
    :param func: If set, filters by the node's function from :mod:`pybel.constants` like
     :data:`pybel.constants.GENE`, :data:`pybel.constants.RNA`,  :data:`pybel.constants.PROTEIN`, or
     :data:`pybel.constants.BIOPROCESS`
    :param prune_threshold: Removes nodes with less than or equal to this number of connections. Defaults to :code:`1`
    :return: An iterable over nodes with only a connection to one node
    """
    for node in graph.nodes(data=True):
        if func and func != node.function:
            continue

        if graph.in_degree(node) + graph.out_degree(node) <= prune_threshold:
            yield node
