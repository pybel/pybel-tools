# -*- coding: utf-8 -*-

from .node_filters import filter_nodes, function_inclusion_filter_builder, namespace_inclusion_builder, \
    function_namespace_inclusion_builder
from .. import pipeline

__all__ = [
    'get_nodes_by_function',
    'get_nodes_by_namespace',
    'get_nodes_by_function_namespace',
]


@pipeline.in_place_mutator
def get_nodes_by_function(graph, function):
    """Get all nodes of a given type.

    :param pybel.BELGraph graph: A BEL graph
    :param str function: The BEL function to filter by
    :return: An iterable of all BEL nodes with the given function
    :rtype: iter[tuple]
    """
    return filter_nodes(graph, function_inclusion_filter_builder(function))


@pipeline.in_place_mutator
def get_nodes_by_namespace(graph, namespace):
    """Returns an iterator over nodes with the given namespace

    :param pybel.BELGraph graph: A BEL graph
    :param str or list[str] namespace: The namespace or list of namespaces to filter
    :return: An iterable over BEL nodes with the given function and namespace
    :rtype: iter[tuple]
    """
    return filter_nodes(graph, namespace_inclusion_builder(namespace))


@pipeline.in_place_mutator
def get_nodes_by_function_namespace(graph, function, namespace):
    """Returns an iterator over nodes with the given function and namespace

    :param pybel.BELGraph graph: A BEL graph
    :param str function: The function to filter
    :param str namespace: The namespace to filter
    :return: An iterable over BEL nodes with the given function and namespace
    :rtype: iter[tuple]
    """
    return filter_nodes(graph, function_namespace_inclusion_builder(function, namespace))
