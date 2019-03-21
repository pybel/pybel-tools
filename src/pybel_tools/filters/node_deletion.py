# -*- coding: utf-8 -*-

"""Node deletion functions to supplement :mod:`pybel.struct.filters`."""

from pybel import BELGraph
from pybel.struct.filters.node_predicate_builders import function_inclusion_filter_builder
from pybel.struct.mutation import remove_filtered_nodes
from pybel.struct.pipeline import in_place_transformation
from pybel.struct.pipeline.decorators import register_deprecated
from pybel.typing import Strings
from .node_filters import function_namespace_inclusion_builder, namespace_inclusion_builder

__all__ = [
    'remove_nodes_by_function',
    'remove_nodes_by_namespace',
    'remove_nodes_by_function_namespace',
]


@in_place_transformation
def remove_nodes_by_function(graph: BELGraph, func: Strings) -> None:
    """Remove nodes with the given function.

    This could be useful directly to remove pathologies.
    """
    remove_filtered_nodes(graph, function_inclusion_filter_builder(func))


@in_place_transformation
def remove_nodes_by_namespace(graph: BELGraph, namespace: Strings) -> None:
    """Remove nodes with the given namespace.

    This might be useful to exclude information learned about distant species, such as excluding all information
    from MGI and RGD in diseases where mice and rats don't give much insight to the human disease mechanism.
    """
    remove_filtered_nodes(graph, namespace_inclusion_builder(namespace))


@register_deprecated('remove_mgi_nodes')
@in_place_transformation
def remove_mouse_nodes(graph: BELGraph) -> None:
    """Remove nodes using the MGI and MGIID namespaces."""
    remove_nodes_by_namespace(graph, ['MGI', 'MGIID'])


@register_deprecated('remove_rgd_nodes')
@in_place_transformation
def remove_rat_nodes(graph: BELGraph) -> None:
    """Remove nodes using the RGD and RGDID namespaces."""
    remove_nodes_by_namespace(graph, ['RGD', 'RGDID'])


@in_place_transformation
def remove_nodes_by_function_namespace(graph: BELGraph, func: str, namespace: Strings) -> None:
    """Remove nodes with the given function and namespace.

    This might be useful to exclude information learned about distant species, such as excluding all information
    from MGI and RGD in diseases where mice and rats don't give much insight to the human disease mechanism.
    """
    remove_filtered_nodes(graph, function_namespace_inclusion_builder(func, namespace))
