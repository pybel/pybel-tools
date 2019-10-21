# -*- coding: utf-8 -*-

"""Node filters to supplement :mod:`pybel.struct.filters.node_filters`."""

from typing import Collection, Iterable, Mapping, Optional, Set, Union

import pybel
from pybel import BELGraph
from pybel.constants import CAUSAL_RELATIONS, FUNCTION, HAS_VARIANT, LABEL, NAMESPACE, PATHOLOGY, RELATION
from pybel.dsl import BaseEntity, Protein, ProteinModification
from pybel.struct.filters import (
    build_node_data_search, build_node_key_search, count_passed_node_filter, data_missing_key_builder,
)
from pybel.struct.filters.typing import NodePredicate, NodePredicates
from pybel.typing import Strings
from ..utils import group_as_sets

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
    'build_node_data_search',
    'build_node_key_search',
    'variants_of',
    'get_variants_to_controllers',
]


def summarize_node_filter(graph: BELGraph, node_filters: NodePredicates) -> None:
    """Print a summary of the number of nodes passing a given set of filters.

    :param graph: A BEL graph
    :param node_filters: A node filter or list/tuple of node filters
    """
    passed = count_passed_node_filter(graph, node_filters)
    print('{}/{} nodes passed'.format(passed, graph.number_of_nodes()))


# Example filters

def node_inclusion_filter_builder(nodes: Iterable[BaseEntity]) -> NodePredicate:
    """Build a filter that only passes on nodes in the given list.

    :param nodes: An iterable of BEL nodes
    """
    node_set = set(nodes)

    def inclusion_filter(_: BELGraph, node: BaseEntity) -> bool:
        """Pass only for a node that is in the enclosed node list.

        :return: If the node is contained within the enclosed node list
        """
        return node in node_set

    return inclusion_filter


def node_exclusion_filter_builder(nodes: Iterable[BaseEntity]) -> NodePredicate:
    """Build a filter that fails on nodes in the given list."""
    node_set = set(nodes)

    def exclusion_filter(_: BELGraph, node: BaseEntity) -> bool:
        """Pass only for a node that isn't in the enclosed node list.

        :return: If the node isn't contained within the enclosed node list
        """
        return node not in node_set

    return exclusion_filter


def _single_function_inclusion_filter_builder(func: str) -> NodePredicate:
    def function_inclusion_filter(_: BELGraph, node: BaseEntity) -> bool:
        """Pass only for a node that has the enclosed function.

        :return: If the node doesn't have the enclosed function
        """
        return node.function == func

    return function_inclusion_filter


def _collection_function_inclusion_builder(funcs: Iterable[str]) -> NodePredicate:
    funcs = set(funcs)

    def functions_inclusion_filter(_: BELGraph, node: BaseEntity) -> bool:
        """Pass only for a node that is one of the enclosed functions.

        :return: If the node doesn't have the enclosed functions
        """
        return node.function in funcs

    return functions_inclusion_filter


def function_inclusion_filter_builder(func: Strings) -> NodePredicate:
    """Build a filter that only passes on nodes of the given function(s).

    :param func: A BEL Function or list/set/tuple of BEL functions
    """
    if isinstance(func, str):
        return _single_function_inclusion_filter_builder(func)

    elif isinstance(func, Iterable):
        return _collection_function_inclusion_builder(func)

    raise ValueError(f'Invalid type for argument: {func}')


def function_exclusion_filter_builder(func: Strings) -> NodePredicate:
    """Build a filter that fails on nodes of the given function(s).

    :param func: A BEL Function or list/set/tuple of BEL functions
    """
    if isinstance(func, str):
        def function_exclusion_filter(_: BELGraph, node: BaseEntity) -> bool:
            """Pass only for a node that doesn't have the enclosed function.

            :return: If the node doesn't have the enclosed function
            """
            return node[FUNCTION] != func

        return function_exclusion_filter

    elif isinstance(func, Iterable):
        functions = set(func)

        def functions_exclusion_filter(_: BELGraph, node: BaseEntity) -> bool:
            """Pass only for a node that doesn't have the enclosed functions.

            :return: If the node doesn't have the enclosed functions
            """
            return node[FUNCTION] not in functions

        return functions_exclusion_filter

    raise ValueError('Invalid type for argument: {}'.format(func))


def function_namespace_inclusion_builder(func: str, namespace: Strings) -> NodePredicate:
    """Build a filter function for matching the given BEL function with the given namespace or namespaces.

    :param func: A BEL function
    :param namespace: The namespace to search by
    """
    if isinstance(namespace, str):
        def function_namespaces_filter(_: BELGraph, node: BaseEntity) -> bool:
            """Pass only for nodes that have the enclosed function and enclosed namespace."""
            if func != node[FUNCTION]:
                return False
            return NAMESPACE in node and node[NAMESPACE] == namespace

    elif isinstance(namespace, Iterable):
        namespaces = set(namespace)

        def function_namespaces_filter(_: BELGraph, node: BaseEntity) -> bool:
            """Pass only for nodes that have the enclosed function and namespace in the enclose set."""
            if func != node[FUNCTION]:
                return False
            return NAMESPACE in node and node[NAMESPACE] in namespaces

    else:
        raise ValueError('Invalid type for argument: {}'.format(namespace))

    return function_namespaces_filter


def data_contains_key_builder(key: str) -> NodePredicate:  # noqa: D202
    """Build a filter that passes only on nodes that have the given key in their data dictionary.

    :param key: A key for the node's data dictionary
    """

    def data_contains_key(_: BELGraph, node: BaseEntity) -> bool:
        """Pass only for a node that contains the enclosed key in its data dictionary.

        :return: If the node contains the enclosed key in its data dictionary
        """
        return key in node

    return data_contains_key


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


def namespace_inclusion_builder(namespace: str) -> NodePredicate:  # noqa: D202
    """Build a function that filters for nods that include a specific namespace."""

    def _has_namespace(_: BELGraph, node: BaseEntity) -> bool:
        return node.get(NAMESPACE) == namespace

    return _has_namespace


def variants_of(
    graph: BELGraph,
    node: Protein,
    modifications: Optional[Set[str]] = None,
) -> Set[Protein]:
    """Return all variants of the given node."""
    if modifications:
        return _get_filtered_variants_of(graph, node, modifications)

    return {
        v
        for u, v, key, data in graph.edges(keys=True, data=True)
        if (
            u == node and
            data[RELATION] == HAS_VARIANT and
            pybel.struct.has_protein_modification(v)
        )
    }


def _get_filtered_variants_of(
    graph: BELGraph,
    node: Protein,
    modifications: Collection[str],
) -> Set[Protein]:
    return {
        v
        for u, v, key, data in graph.edges(keys=True, data=True)
        if (
            u == node and
            data[RELATION] == HAS_VARIANT and
            pybel.struct.has_protein_modification(v) and
            any(
                variant['identifier']['name'] in modifications
                for variant in v.variants if 'identifier' in variant
                if isinstance(variant, ProteinModification)
            )
        )
    }


def get_variants_to_controllers(
    graph: BELGraph,
    node: Protein,
    modifications: Optional[Set[str]] = None,
    relations: Union[None, str, Set[str]] = None,
) -> Mapping[Protein, Set[Protein]]:
    """Get a mapping from variants of the given node to all of its upstream controllers."""
    variants = variants_of(graph, node, modifications)

    if relations is None:
        relations = CAUSAL_RELATIONS
    elif isinstance(relations, str):
        relations = {relations}

    return group_as_sets(
        (variant, controller)
        for controller, variant, data in graph.in_edges(variants, data=True)
        if data[RELATION] in relations
    )
