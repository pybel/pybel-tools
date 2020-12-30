# -*- coding: utf-8 -*-

"""Node filters to supplement :mod:`pybel.struct.filters.node_filters`."""

import warnings
from typing import Collection, Mapping, Optional, Set, Union

import pybel
import pybel.struct.filters
from pybel import BELGraph
from pybel.constants import CAUSAL_RELATIONS, HAS_VARIANT, RELATION
from pybel.dsl import BaseEntity, Protein, ProteinModification
from pybel.struct.filters import (
    NodePredicate, build_node_data_search, build_node_key_search,
    concatenate_node_predicates, data_missing_key_builder,
)
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


def function_inclusion_filter_builder(func: Strings) -> NodePredicate:
    """Build a filter that only passes on nodes of the given function(s).

    :param func: A BEL Function or list/set/tuple of BEL functions
    """
    warnings.warn('use pybel.struct.function_inclusion_filter_builder', DeprecationWarning)
    return pybel.struct.filters.function_inclusion_filter_builder(func)


def function_exclusion_filter_builder(func: Strings) -> NodePredicate:
    """Build a filter that fails on nodes of the given function(s).

    :param func: A BEL Function or list/set/tuple of BEL functions
    """
    warnings.warn('use pybel.struct.function_exclusion_filter_builder', DeprecationWarning)
    return pybel.struct.filters.function_exclusion_filter_builder(func)


def function_namespace_inclusion_builder(func: Strings, namespace: Strings) -> NodePredicate:
    """Build a filter function for matching the given BEL function with the given namespace or namespaces.

    :param func: A BEL function
    :param namespace: The namespace to search by
    """
    return concatenate_node_predicates([
        pybel.struct.filters.function_inclusion_filter_builder(func),
        pybel.struct.filters.namespace_inclusion_builder(namespace),
    ])


def namespace_inclusion_builder(namespace: str) -> NodePredicate:  # noqa: D202
    """Build a function that filters for nods that include a specific namespace."""
    warnings.warn('use pybel.struct.namespace_inclusion_builder', DeprecationWarning)
    return pybel.struct.filters.namespace_inclusion_builder(namespace)


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
