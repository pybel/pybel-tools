# -*- coding: utf-8 -*-

"""
Node Filters
------------

A node filter is a function that takes two arguments: a :class:`pybel.BELGraph` and a node tuple. It returns a boolean
representing whether the node passed the given test.

This module contains a set of default functions for filtering lists of nodes and building node filtering functions.

A general use for a node filter function is to use the built-in :func:`filter` in code like
:code:`filter(your_node_filter, graph)`
"""

from collections import Iterable

from pybel.constants import FUNCTION, FUSION, HAS_MEMBER, LABEL, NAMESPACE, PATHOLOGY, PROTEIN, RELATION, VARIANTS
from pybel.struct.filters import build_node_data_search, build_node_key_search, data_missing_key_builder
from pybel.struct.filters.node_filters import count_passed_node_filter, filter_nodes
from ..constants import CNAME

__all__ = [
    'summarize_node_filter',
    'node_inclusion_filter_builder',
    'node_exclusion_filter_builder',
    'function_inclusion_filter_builder',
    'function_exclusion_filter_builder',
    'function_namespace_inclusion_builder',
    'namespace_inclusion_builder',
    'data_contains_key_builder',
    'data_missing_key_builder',
    'node_has_cname',
    'node_missing_cname',
    'node_has_label',
    'node_missing_label',
    'include_pathology_filter',
    'exclude_pathology_filter',
    'build_node_data_search',
    'build_node_key_search',
    'build_node_cname_search',
    'iter_undefined_families',
]


def summarize_node_filter(graph, node_filters):
    """Prints a summary of the number of nodes passing a given set of filters

    :param pybel.BELGraph graph: A BEL graph
    :param node_filters: A node filter or list/tuple of node filters
    :type node_filters: types.FunctionType or iter[types.FunctionType]
    """
    passed = count_passed_node_filter(graph, node_filters)
    print('{}/{} nodes passed'.format(passed, graph.number_of_nodes()))


# Example filters

def node_inclusion_filter_builder(nbunch):
    """Builds a filter that only passes on nodes in the given list

    :param iter[tuple] nbunch: An iterable of BEL nodes
    :return: A node filter (graph, node) -> bool
    :rtype: types.FunctionType
    """
    node_set = set(nbunch)

    def inclusion_filter(graph, node):
        """Passes only for a node that is in the enclosed node list

        :param pybel.BELGraph graph: A BEL Graph
        :param tuple node: A BEL node
        :return: If the node is contained within the enclosed node list
        :rtype: bool
        """
        return node in node_set

    return inclusion_filter


def node_exclusion_filter_builder(nodes):
    """Builds a filter that fails on nodes in the given list

    :param nodes: A list of nodes
    :type nodes: list
    :return: A node filter (graph, node) -> bool
    :rtype: types.FunctionType
    """
    node_set = set(nodes)

    def exclusion_filter(graph, node):
        """Passes only for a node that isn't in the enclosed node list

        :param pybel.BELGraph graph: A BEL Graph
        :param tuple node: A BEL node
        :return: If the node isn't contained within the enclosed node list
        :rtype: bool
        """
        return node not in node_set

    return exclusion_filter


def _single_function_inclusion_filter_builder(func):
    def function_inclusion_filter(graph, node):
        """Pass only for a node that has the enclosed function.

        :param BELGraph graph: A BEL Graph
        :param tuple node: A BEL node
        :return: If the node doesn't have the enclosed function
        :rtype: bool
        """
        return graph.node[node][FUNCTION] == func

    return function_inclusion_filter


def _collection_function_inclusion_builder(funcs):
    funcs = set(funcs)

    def functions_inclusion_filter(graph, node):
        """Pass only for a node that is one of the enclosed functions.

        :param BELGraph graph: A BEL Graph
        :param tuple node: A BEL node
        :return: If the node doesn't have the enclosed functions
        :rtype: bool
        """
        return graph.node[node][FUNCTION] in funcs

    return functions_inclusion_filter


def function_inclusion_filter_builder(func):
    """Build a filter that only passes on nodes of the given function(s).

    :param func: A BEL Function or list/set/tuple of BEL functions
    :type func: str or iter[str]
    :return: A node predicate
    :rtype: (pybel.BELGraph, tuple) -> bool
    """
    if isinstance(func, str):
        return _single_function_inclusion_filter_builder(func)

    elif isinstance(func, Iterable):
        return _collection_function_inclusion_builder(func)

    raise ValueError('Invalid type for argument: {}'.format(func))


def function_exclusion_filter_builder(func):
    """Build a filter that fails on nodes of the given function(s).

    :param func: A BEL Function or list/set/tuple of BEL functions
    :type func: str or list[str] or tuple[str] or set[str]
    :return: A node predicate
    :rtype: (pybel.BELGraph, tuple) -> bool
    """
    if isinstance(func, str):
        def function_exclusion_filter(graph, node):
            """Pass only for a node that doesn't have the enclosed function.

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple node: A BEL node
            :return: If the node doesn't have the enclosed function
            :rtype: bool
            """
            return graph.node[node][FUNCTION] != func

        return function_exclusion_filter

    elif isinstance(func, Iterable):
        functions = set(func)

        def functions_exclusion_filter(graph, node):
            """Pass only for a node that doesn't have the enclosed functions.

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple node: A BEL node
            :return: If the node doesn't have the enclosed functions
            :rtype: bool
            """
            return graph.node[node][FUNCTION] not in functions

        return functions_exclusion_filter

    raise ValueError('Invalid type for argument: {}'.format(func))


def function_namespace_inclusion_builder(func, namespace):
    """Build a filter function for matching the given BEL function with the given namespace or namespaces.

    :param str func: A BEL function
    :param str or iter[str] namespace: The namespace to serach by
    :return: A node predicate
    :rtype: (pybel.BELGraph, tuple) -> bool
    """
    if isinstance(namespace, str):
        def function_namespace_filter(graph, node):
            """Passes only for nodes that have the enclosed function and enclosed namespace

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple node: A BEL node
            :rtype: bool
            """
            if func != graph.node[node][FUNCTION]:
                return False
            return NAMESPACE in graph.node[node] and graph.node[node][NAMESPACE] == namespace

        return function_namespace_filter

    elif isinstance(namespace, Iterable):
        namespaces = set(namespace)

        def function_namespaces_filter(graph, node):
            """Passes only for nodes that have the enclosed function and namespace in the enclose set

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple node: A BEL node
            :rtype: bool
            """
            if func != graph.node[node][FUNCTION]:
                return False
            return NAMESPACE in graph.node[node] and graph.node[node][NAMESPACE] in namespaces

        return function_namespaces_filter

    raise ValueError('Invalid type for argument: {}'.format(namespace))


def namespace_inclusion_builder(namespace):
    """Build a predicate for namespace inclusion.

    :param str or iter[str] namespace: A namespace or iter of namespaces
    :return: A node predicate
    :rtype: (pybel.BELGraph, tuple) -> bool
    """
    if isinstance(namespace, str):

        def namespace_filter(graph, node):
            """Passes only for a node that has the enclosed namespace

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple node: A BEL node
            :rtype: bool
            """
            return NAMESPACE in graph.node[node] and graph.node[node][NAMESPACE] == namespace

        return namespace_filter

    elif isinstance(namespace, Iterable):
        namespaces = set(namespace)

        def namespaces_filter(graph, node):
            """Pass only for a node that has a namespace in the enclosed set.

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple node: A BEL node
            :rtype: bool
            """
            return NAMESPACE in graph.node[node] and graph.node[node][NAMESPACE] in namespaces

        return namespaces_filter

    raise ValueError('Invalid type for argument: {}'.format(namespace))


def data_contains_key_builder(key):
    """Build a filter that passes only on nodes that have the given key in their data dictionary.

    :param str key: A key for the node's data dictionary
    :return: A node predicate
    :rtype: (pybel.BELGraph, tuple) -> bool
    """

    def data_contains_key(graph, node):
        """Pass only for a node that contains the enclosed key in its data dictionary.

        :param pybel.BELGraph graph: A BEL Graph
        :param tuple node: A BEL node
        :return: If the node contains the enclosed key in its data dictionary
        :rtype: bool
        """
        return key in graph.node[node]

    return data_contains_key


#: Passes for nodes that have been annotated with a canonical name
node_has_cname = data_contains_key_builder(CNAME)

#: Fails for nodes that have been annotated with a canonical name
node_missing_cname = data_missing_key_builder(CNAME)

#: Passes for nodes that have been annotated with a label
node_has_label = data_contains_key_builder(LABEL)

#: Fails for nodes that have been annotated with a label
node_missing_label = data_missing_key_builder(LABEL)

# Default Filters

#: A filter that passes for nodes that are :data:`pybel.constants.PATHOLOGY`
include_pathology_filter = function_inclusion_filter_builder(PATHOLOGY)

#: A filter that fails for nodes that are :data:`pybel.constants.PATHOLOGY`
exclude_pathology_filter = function_exclusion_filter_builder(PATHOLOGY)


# TODO node filter that is false for abundances with no in-edges


def build_node_cname_search(query):
    """Search nodes' canonical names. Is a thin wrapper around :func:`build_node_key_search`.

    :param query: The query string or strings to check if they're in the node name
    :type query: str or iter[str]
    :return: A node predicate
    :rtype: (pybel.BELGraph, tuple) -> bool
    """
    return build_node_key_search(query, CNAME)


def iter_undefined_families(graph, namespace):
    """Finds protein families from a given namespace (Such as SFAM) that aren't qualified by members

    :param pybel.BELGraph graph: A BEL graph
    :param namespace: The namespace to filter by
    :type namespace: str or iter[str]
    :return: An iterator over nodes that don't
    :rtype: iter[tuple]
    """
    node_predicates = [
        function_inclusion_filter_builder(PROTEIN),
        namespace_inclusion_builder(namespace)
    ]

    for node in filter_nodes(graph, node_predicates):
        if VARIANTS or FUSION in graph.node[node]:
            continue

        relations = {
            d[RELATION]
            for _, v, d in graph.out_edges_iter(node, data=True)
        }

        if HAS_MEMBER in relations:
            continue

        yield node
