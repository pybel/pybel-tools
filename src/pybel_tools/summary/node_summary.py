# -*- coding: utf-8 -*-

"""This module contains functions that provide summaries of the nodes in a graph"""

from collections import Counter, defaultdict

from pybel.constants import *

__all__ = [
    'count_functions',
    'get_functions',
    'count_namespaces',
    'get_namespaces',
    'count_names',
    'get_names_by_namespace',
    'get_unused_namespaces',
    'get_names',
]


def count_functions(graph):
    """Counts the frequency of each function present in a graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A Counter from {function: frequency}
    :rtype: collections.Counter
    """
    return Counter(
        data[FUNCTION]
        for _, data in graph.nodes_iter(data=True)
    )


def get_functions(graph):
    """Gets the set of all functions used in this graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A set of functions
    :rtype: set[str]
    """
    return set(count_functions(graph))


def count_namespaces(graph):
    """Counts the frequency of each namespace across all nodes (that have namespaces)

    :param pybel.BELGraph graph: A BEL graph
    :return: A Counter from {namespace: frequency}
    :rtype: collections.Counter
    """
    return Counter(
        data[NAMESPACE]
        for _, data in graph.nodes_iter(data=True)
        if NAMESPACE in data
    )


def get_namespaces(graph):
    """Gets the set of all namespaces used in this graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A set of namespaces
    :rtype: set[str]
    """
    return set(count_namespaces(graph))


def get_unused_namespaces(graph):
    """Gets the set of all namespaces that are defined in a graph, but are never used.
    
    :param pybel.BELGraph graph: A BEL graph
    :return: A set of namespaces that are included but not used
    :rtype: set[str] 
    """
    defined_namespaces = set(graph.namespace_pattern) | set(graph.namespace_url) | set(graph.namespace_owl)
    return defined_namespaces - get_namespaces(graph)


def count_names(graph):
    """Counts all names through the graph by the NAME tag in the nodes' data dictionaries.

    This is useful to identify which nodes appear with the same name in multiple namespaces, or to identify variants

    :param pybel.BELGraph graph: A BEL graph
    :return: A Counter from {names: frequency}
    :rtype: collections.Counter
    """
    return Counter(
        data[NAME]
        for _, data in graph.nodes_iter(data=True)
        if NAME in data
    )


def get_names_by_namespace(graph, namespace):
    """Get the set of all of the names in a given namespace that are in the graph

    :param pybel.BELGraph graph: A BEL graph
    :param str namespace: A namespace
    :return: A set of names belonging to the given namespace that are in the given graph
    :rtype: set[str]
    """
    result = set()

    for _, data in graph.nodes_iter(data=True):
        if NAMESPACE in data and data[NAMESPACE] == namespace:
            result.add(data[NAME])
        elif FUSION in data:
            if data[FUSION][PARTNER_3P][NAMESPACE] == namespace:
                result.add(data[FUSION][PARTNER_3P][NAME])

            if data[FUSION][PARTNER_5P][NAMESPACE] == namespace:
                result.add(data[FUSION][PARTNER_5P][NAME])

    return result


def get_names(graph):
    """Get a dictionary of {namespace: set of names} present in the graph.

    Roughly equivalent to:

    >>> {namespace: get_names_by_namespace(graph, namespace) for namespace in get_namespaces(graph)}

    :param pybel.BELGraph graph: A BEL graph
    :return: A dictionary of {namespace: set of names}
    :rtype: dict[str, set[str]]
    """
    result = defaultdict(set)

    for _, data in graph.nodes_iter(data=True):
        if NAMESPACE in data:
            result[data[NAMESPACE]].add(data[NAME])
        elif FUSION in data:
            result[data[FUSION][PARTNER_3P][NAMESPACE]].add(data[FUSION][PARTNER_3P][NAME])
            result[data[FUSION][PARTNER_5P][NAMESPACE]].add(data[FUSION][PARTNER_5P][NAME])

    return dict(result)
