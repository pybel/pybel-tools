# -*- coding: utf-8 -*-

from pybel.struct.filters import filter_nodes
from .node_filters import (
    function_inclusion_filter_builder,
    namespace_inclusion_builder,
    function_namespace_inclusion_builder
)
from .. import pipeline

__all__ = [
    'remove_filtered_nodes',
    'remove_nodes_by_function',
    'remove_nodes_by_namespace',
    'remove_nodes_by_function_namespace',
]


@pipeline.in_place_mutator
def remove_filtered_nodes(graph, node_filters):
    """Removes nodes passing the given node filters

    :param pybel.BELGraph graph: A BEL graph
    :param node_filters: A node filter list of node filters (graph, node) -> bool
    :type node_filters: types.FunctionType or iter[types.FunctionType]
    """
    nodes = list(filter_nodes(graph, node_filters))
    graph.remove_nodes_from(nodes)


@pipeline.in_place_mutator
def remove_nodes_by_function(graph, function):
    """Removes nodes with the given function.

    This could be useful directly to remove pathologies.
    
    :param pybel.BELGraph graph: A BEL graph
    :param str function: The function to filter
    """
    remove_filtered_nodes(graph, function_inclusion_filter_builder(function))


@pipeline.in_place_mutator
def remove_nodes_by_namespace(graph, namespace):
    """Removes nodes with the given  namespace.

    This might be useful to exclude information learned about distant species, such as excluding all information
    from MGI and RGD in diseases where mice and rats don't give much insight to the human disease mechanism.

    :param pybel.BELGraph graph: A BEL graph
    :param str or iter[str] namespace: The namespace to filter or iterable of namespaces
    """
    remove_filtered_nodes(graph, namespace_inclusion_builder(namespace))


@pipeline.in_place_mutator
def remove_mgi_nodes(graph):
    """Removes MGI nodes.

    :param pybel.BELGraph graph: A BEL graph
    """
    remove_nodes_by_namespace(graph, 'MGI')


@pipeline.in_place_mutator
def remove_rgd_nodes(graph):
    """Removes RGD nodes.

    :param pybel.BELGraph graph: A BEL graph
    """
    remove_nodes_by_namespace(graph, 'RGD')


@pipeline.in_place_mutator
def remove_nodes_by_function_namespace(graph, function, namespace):
    """Removes nodes with the given function and namespace.

    This might be useful to exclude information learned about distant species, such as excluding all information
    from MGI and RGD in diseases where mice and rats don't give much insight to the human disease mechanism.

    :param pybel.BELGraph graph: A BEL graph
    :param str function: The function to filter
    :param str or iter[str] namespace: The namespace to filter or iterable of namespaces
    """
    remove_filtered_nodes(graph, function_namespace_inclusion_builder(function, namespace))
