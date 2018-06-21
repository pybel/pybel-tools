# -*- coding: utf-8 -*-

from pybel.struct import filter_nodes, get_nodes_by_function
from .node_filters import function_namespace_inclusion_builder, namespace_inclusion_builder

__all__ = [
    'get_nodes_by_function',
    'get_nodes_by_namespace',
    'get_nodes_by_function_namespace',
]


def get_nodes_by_namespace(graph, namespace):
    """Returns an iterator over nodes with the given namespace

    :param pybel.BELGraph graph: A BEL graph
    :param str or list[str] namespace: The namespace or list of namespaces to filter
    :return: An iterable over BEL nodes with the given function and namespace
    :rtype: iter[tuple]
    """
    return filter_nodes(graph, namespace_inclusion_builder(namespace))


def get_nodes_by_function_namespace(graph, func, namespace):
    """Returns an iterator over nodes with the given function and namespace

    :param pybel.BELGraph graph: A BEL graph
    :param str func: The function to filter
    :param str or iter[str] namespace: The namespace to filter
    :return: An iterable over BEL nodes with the given function and namespace
    :rtype: iter[tuple]
    """
    return filter_nodes(graph, function_namespace_inclusion_builder(func, namespace))
