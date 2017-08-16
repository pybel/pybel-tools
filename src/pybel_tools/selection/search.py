# -*- coding: utf-8 -*-

from ..filters.node_filters import filter_nodes, build_node_name_search, build_node_cname_search

__all__ = [
    'search_node_names',
    'search_node_cnames'
]


def search_node_names(graph, query):
    """Searches for nodes containing a given string

    :param pybel.BELGraph graph: A BEL graph
    :param query: The search query
    :type query: str or iter[str]
    :return: An iterator over nodes whose names match the search query
    :rtype: iter
    """

    return filter_nodes(graph, build_node_name_search(query))


def search_node_cnames(graph, query):
    """Searches for nodes whose canonical names contain a given string(s)

    :param pybel.BELGraph graph: A BEL graph
    :param query: The search query
    :type query: str or iter[str]
    :return: An iterator over nodes whose canonical names match the search query
    :rtype: iter
    """
    return filter_nodes(graph, build_node_cname_search(query))
