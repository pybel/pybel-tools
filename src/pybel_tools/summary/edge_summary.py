# -*- coding: utf-8 -*-

"""This module contains functions that provide summaries of the edges in a graph."""

import itertools as itt
from collections import Counter, defaultdict
from typing import Iterable, List, Mapping, Optional, Set, Tuple, TypeVar

from pybel import BELGraph
from pybel.constants import ANNOTATIONS, RELATION
from pybel.dsl import BaseEntity
from pybel.struct.filters.edge_predicates import edge_has_annotation
from pybel.struct.filters.typing import NodePredicate
from pybel.struct.summary import (
    count_annotations, count_pathologies, count_relations, get_annotations, get_unused_annotations,
    get_unused_list_annotation_values, iter_annotation_value_pairs, iter_annotation_values,
)
from .contradictions import pair_has_contradiction

__all__ = [
    'count_relations',
    'get_edge_relations',
    'count_unique_relations',
    'count_annotations',
    'get_annotations',
    'get_annotations_containing_keyword',
    'count_annotation_values',
    'count_annotation_values_filtered',
    'pair_is_consistent',
    'get_consistent_edges',
    'get_contradictory_pairs',
    'count_pathologies',
    'get_unused_annotations',
    'get_unused_list_annotation_values',
]

A = TypeVar('A')
B = TypeVar('B')


def group_dict_set(iterator: Iterable[Tuple[A, B]]) -> Mapping[A, Set[B]]:
    """Make a dict that accumulates the values for each key in an iterator of doubles."""
    d = defaultdict(set)
    for key, value in iterator:
        d[key].add(value)
    return dict(d)


def get_edge_relations(graph: BELGraph) -> Mapping[Tuple[BaseEntity, BaseEntity], Set[str]]:
    """Build a dictionary of {node pair: set of edge types}."""
    return group_dict_set(
        ((u, v), d[RELATION])
        for u, v, d in graph.edges(data=True)
    )


def count_unique_relations(graph: BELGraph) -> Counter:
    """Return a histogram of the different types of relations present in a graph.

    Note: this operation only counts each type of edge once for each pair of nodes
    """
    return Counter(itt.chain.from_iterable(get_edge_relations(graph).values()))


def get_annotations_containing_keyword(graph: BELGraph, keyword: str) -> List[Mapping[str, str]]:
    """Get annotation/value pairs for values for whom the search string is a substring.

    :param graph: A BEL graph
    :param keyword: Search for annotations whose values have this as a substring
    """
    return [
        {
            'annotation': annotation,
            'value': value
        }
        for annotation, value in iter_annotation_value_pairs(graph)
        if keyword.lower() in value.lower()
    ]


def count_annotation_values(graph: BELGraph, annotation: str) -> Counter:
    """Count in how many edges each annotation appears in a graph.

    :param graph: A BEL graph
    :param annotation: The annotation to count
    :return: A Counter from {annotation value: frequency}
    """
    return Counter(iter_annotation_values(graph, annotation))


def count_annotation_values_filtered(
        graph: BELGraph,
        annotation: str,
        source_predicate: Optional[NodePredicate] = None,
        target_predicate: Optional[NodePredicate] = None,
) -> Counter:
    """Count in how many edges each annotation appears in a graph, but filter out source nodes and target nodes.

    See :func:`pybel_tools.utils.keep_node` for a basic filter.

    :param graph: A BEL graph
    :param annotation: The annotation to count
    :param source_predicate: A predicate (graph, node) -> bool for keeping source nodes
    :param target_predicate: A predicate (graph, node) -> bool for keeping target nodes
    :return: A Counter from {annotation value: frequency}
    """
    if source_predicate and target_predicate:
        return Counter(
            data[ANNOTATIONS][annotation]
            for u, v, data in graph.edges(data=True)
            if edge_has_annotation(data, annotation) and source_predicate(graph, u) and target_predicate(graph, v)
        )
    elif source_predicate:
        return Counter(
            data[ANNOTATIONS][annotation]
            for u, v, data in graph.edges(data=True)
            if edge_has_annotation(data, annotation) and source_predicate(graph, u)
        )
    elif target_predicate:
        return Counter(
            data[ANNOTATIONS][annotation]
            for u, v, data in graph.edges(data=True)
            if edge_has_annotation(data, annotation) and target_predicate(graph, u)
        )
    else:
        return Counter(
            data[ANNOTATIONS][annotation]
            for u, v, data in graph.edges(data=True)
            if edge_has_annotation(data, annotation)
        )


def pair_is_consistent(graph: BELGraph, u: BaseEntity, v: BaseEntity) -> Optional[str]:
    """Return if the edges between the given nodes are consistent, meaning they all have the same relation.

    :return: If the edges aren't consistent, return false, otherwise return the relation type
    """
    relations = {data[RELATION] for data in graph[u][v].values()}

    if 1 != len(relations):
        return

    return list(relations)[0]


def get_contradictory_pairs(graph: BELGraph) -> Iterable[Tuple[BaseEntity, BaseEntity]]:
    """Iterate over contradictory node pairs in the graph based on their causal relationships.

    :return: An iterator over (source, target) node pairs that have contradictory causal edges
    """
    for u, v in graph.edges():
        if pair_has_contradiction(graph, u, v):
            yield u, v


def get_consistent_edges(graph: BELGraph) -> Iterable[Tuple[BaseEntity, BaseEntity]]:
    """Yield pairs of (source node, target node) for which all of their edges have the same type of relation.

    :return: An iterator over (source, target) node pairs corresponding to edges with many inconsistent relations
    """
    for u, v in graph.edges():
        if pair_is_consistent(graph, u, v):
            yield u, v
