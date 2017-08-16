# -*- coding: utf-8 -*-

"""

A metapath can be defined with two levels of granularity:

- Low: A list of BEL functions representing the types of entities in a given path
- High: An alternating list of BEL functions and BEL relations representing the types of entities in a given path and
  their relations

"""

from functools import lru_cache

from pybel.constants import FUNCTION

__all__ = [
    'convert_path_to_metapath',
    'get_walks_exhaustive',
    'match_simple_metapath',
]


def convert_path_to_metapath(graph, nodes):
    """Converts a list of nodes to their corresponding functions

    :param list[tuple] nodes: A list of BEL node tuples
    :rtype: list[str]
    """
    return [
        graph.node[node][FUNCTION]
        for node in nodes
    ]


@lru_cache(maxsize=None)
def get_walks_exhaustive(graph, node, length):
    """Gets all walks under a given length starting at a given node

    :param networkx.Graph graph: A graph
    :param node: Starting node
    :param int length: The length of walks to get
    :return: A list of paths
    :rtype: list[tuple]
    """
    if 0 == length:
        return (node,),

    return tuple(
        (node, key) + path
        for neighbor in graph.edge[node]
        for path in get_walks_exhaustive(graph, neighbor, length - 1)
        if node not in path
        for key in graph.edge[node][neighbor]
    )


def match_simple_metapath(graph, node, simple_metapath):
    """Matches a simple metapath starting at the given node

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :param list[str] simple_metapath: A list of BEL Functions
    :return: An iterable over paths from the node matching the metapath
    :rtype: iter[tuple]
    """
    if 0 == len(simple_metapath):
        yield node,

    else:
        for neighbor in graph.edge[node]:
            if graph.node[neighbor][FUNCTION] == simple_metapath[0]:
                for path in match_simple_metapath(graph, neighbor, simple_metapath[1:]):
                    if node not in path:
                        yield (node,) + path


def convert_simple_walk(graph, simple_walk):
    """Converts a walk into a sequence of BEL functions

    :param pybel.BELGraph graph: A BEL graph
    :param iter[tuple] simple_walk: An iterable of BEL nodes
    :return: A list of BEL functions of the walk
    :rtype: list[str]
    """
    return [
        graph.node[node][FUNCTION]
        for node in simple_walk
    ]


def match_complex_metapath(graph, node, complex_metapath):
    """Matches a complex metapath starting at the given node

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :param list[str] complex_metapath: An iterable of alternating BEL nodes and relations
    :return: An iterable over paths from the node matching the metapath
    :rtype: iter[tuple]
    """
    raise NotImplemented


def convert_complex_walk(graph, complex_walk):
    """Converts a walk into an alternative sequence of BEL functions and BEL relations, starting and ending
    with a BEL function

    :param pybel.BELGraph graph: A BEL graph
    :param iter[tuple] complex_walk: An iterable of alternating BEL nodes and relations
    :return: An alternating list of BEL functions and relations of the walk
    :rtype: list[str]
    """
    raise NotImplemented
