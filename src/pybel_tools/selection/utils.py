# -*- coding: utf-8 -*-

from pybel.constants import FUNCTION

__all__ = [
    'get_leaves_by_type',
]


def get_leaves_by_type(graph, function=None, prune_threshold=1):
    """Returns an iterable over all nodes in graph (in-place) with only a connection to one node. Useful for gene and
     RNA. Allows for optional filter by function type.

    :param pybel.BELGraph graph: A BEL graph
    :param function: If set, filters by the node's function from :mod:`pybel.constants` like 
                    :data:`pybel.constants.GENE`, :data:`pybel.constants.RNA`,  :data:`pybel.constants.PROTEIN`, or 
                    :data:`pybel.constants.BIOPROCESS`
    :type function: str
    :param prune_threshold: Removes nodes with less than or equal to this number of connections. Defaults to :code:`1`
    :type prune_threshold: int
    :return: An iterable over nodes with only a connection to one node
    :rtype: iter
    """
    for node, data in graph.nodes_iter(data=True):
        if function and function != data.get(FUNCTION):
            continue

        if graph.in_degree(node) + graph.out_degree(node) <= prune_threshold:
            yield node
