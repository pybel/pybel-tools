# -*- coding: utf-8 -*-

"""Node filters to supplement :mod:`pybel.struct.filters.node_filters`."""

from typing import Collection, Iterable, Mapping, Optional, Set, Union

import pybel
from pybel import BELGraph, BaseAbundance
from pybel.constants import CAUSAL_RELATIONS, HAS_VARIANT, RELATION
from pybel.dsl import BaseEntity, Protein, ProteinModification
from pybel.struct.filters import build_node_data_search, build_node_key_search, data_missing_key_builder
from pybel.struct.filters.typing import NodePredicate
from pybel.typing import Strings
from ..utils import group_as_sets

__all__ = [
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
            return node.function != func

        return function_exclusion_filter

    elif isinstance(func, Iterable):
        functions = set(func)

        def functions_exclusion_filter(_: BELGraph, node: BaseEntity) -> bool:
            """Pass only for a node that doesn't have the enclosed functions.

            :return: If the node doesn't have the enclosed functions
            """
            return node.function not in functions

        return functions_exclusion_filter

    raise ValueError('Invalid type for argument: {}'.format(func))


def function_namespace_inclusion_builder(func: str, namespace: Strings) -> NodePredicate:
    """Build a filter function for matching the given BEL function with the given namespace or namespaces.

    :param func: A BEL function
    :param namespace: The namespace to search by
    """
    if isinstance(namespace, str):
        namespace = namespace.lower()

        def function_namespaces_filter(_: BELGraph, node: BaseEntity) -> bool:
            """Pass only for nodes that have the enclosed function and enclosed namespace."""
            return (
                node.function.lower() == func.lower()
                and isinstance(node, BaseAbundance)
                and node.namespace == namespace
            )

    elif isinstance(namespace, Iterable):
        namespaces = {n.lower() for n in namespace}

        def function_namespaces_filter(_: BELGraph, node: BaseEntity) -> bool:
            """Pass only for nodes that have the enclosed function and namespace in the enclose set."""
            return (
                node.function == func
                and isinstance(node, BaseAbundance)
                and node.namespace.lower() in namespaces
            )

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


def namespace_inclusion_builder(namespace: str) -> NodePredicate:  # noqa: D202
    """Build a function that filters for nods that include a specific namespace."""

    def _has_namespace(_: BELGraph, node: BaseEntity) -> bool:
        return isinstance(node, BaseAbundance) and node.namespace == namespace

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
        target
        for source, target, key, data in graph.edges(keys=True, data=True)
        if (
            source == node
            and data[RELATION] == HAS_VARIANT
            and pybel.struct.has_protein_modification(target)
            and any(
                variant.name in modifications
                for variant in target.variants
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
