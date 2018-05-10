# -*- coding: utf-8 -*-

from pybel.struct.filters import filter_nodes
from ..filters.node_filters import build_node_cname_search, build_node_name_search

__all__ = [
    'search_node_names',
    'search_node_cnames'
]


def search_node_names(graph, query):
    """Search for nodes containing a given string(s).

    :param pybel.BELGraph graph: A BEL graph
    :param query: The search query
    :type query: str or iter[str]
    :return: An iterator over nodes whose names match the search query
    :rtype: iter

    Example:

    .. code-block:: python

        >>> from pybel.examples import sialic_acid_graph
        >>> from pybel_tools.selection import search_node_names
        >>> list(search_node_names(sialic_acid_graph, 'CD33'))
        [('Protein', 'HGNC', 'CD33'), ('Protein', 'HGNC', 'CD33', ('pmod', ('bel', 'Ph')))]
    """
    return filter_nodes(graph, build_node_name_search(query))


def search_node_cnames(graph, query):
    """Search for nodes whose canonical names contain a given string(s).

    :param pybel.BELGraph graph: A BEL graph
    :param query: The search query
    :type query: str or iter[str]
    :return: An iterator over nodes whose canonical names match the search query
    :rtype: iter
    """
    return filter_nodes(graph, build_node_cname_search(query))
