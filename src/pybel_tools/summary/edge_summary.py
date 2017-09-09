# -*- coding: utf-8 -*-

"""This module contains functions that provide summaries of the edges in a graph"""

import itertools as itt

from collections import Counter, defaultdict

from pybel.constants import (
    RELATION,
    ANNOTATIONS,
    FUNCTION,
    PATHOLOGY,
    CAUSAL_INCREASE_RELATIONS,
    CAUSAL_DECREASE_RELATIONS,
    CAUSES_NO_CHANGE,
)
from pybel.struct.filters.node_filters import keep_node_permissive
from ..utils import get_value_sets, check_has_annotation, graph_edge_data_iter

__all__ = [
    'count_relations',
    'get_edge_relations',
    'count_unique_relations',
    'count_annotations',
    'get_annotations',
    'get_annotation_values_by_annotation',
    'get_annotations_containing_keyword',
    'count_annotation_values',
    'get_annotation_values',
    'count_annotation_values_filtered',
    'get_all_relations',
    'pair_is_consistent',
    'get_consistent_edges',
    'pair_has_contradiction',
    'get_inconsistent_edges',
    'get_contradictory_pairs',
    'count_pathologies',
    'get_tree_annotations',
    'relation_set_has_contradictions',
    'get_unused_annotations',
    'get_unused_list_annotation_values',
]


def count_relations(graph):
    """Returns a histogram over all relationships in a graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A Counter from {relation type: frequency}
    :rtype: collections.Counter
    """
    return Counter(d[RELATION] for _, _, d in graph.edges_iter(data=True))


def get_edge_relations(graph):
    """Builds a dictionary of {node pair: set of edge types}

    :param pybel.BELGraph graph: A BEL graph
    :return: A dictionary of {(node, node): set of edge types}
    :rtype: dict
    """
    edge_relations = defaultdict(set)
    for u, v, d in graph.edges_iter(data=True):
        edge_relations[u, v].add(d[RELATION])
    return dict(edge_relations)


def count_unique_relations(graph):
    """Returns a histogram of the different types of relations present in a graph.

    Note: this operation only counts each type of edge once for each pair of nodes

    :param pybel.BELGraph graph: A BEL graph
    :return: Counter from {relation type: frequency}
    :rtype: collections.Counter
    """
    return Counter(itt.chain.from_iterable(get_edge_relations(graph).values()))


def count_annotations(graph):
    """Counts how many times each annotation is used in the graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A Counter from {annotation key: frequency}
    :rtype: collections.Counter
    """
    return Counter(
        key
        for data in graph_edge_data_iter(graph)
        if ANNOTATIONS in data
        for key in data[ANNOTATIONS]
    )


def get_annotations(graph):
    """Gets the set of annotations used in the graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A set of annotation keys
    :rtype: set[str]
    """
    return set(count_annotations(graph))


def get_unused_annotations(graph):
    """Gets the set of all annotations that are defined in a graph, but are never used.

    :param pybel.BELGraph graph: A BEL graph
    :return: A set of annotations
    :rtype: set[str] 
    """
    defined_annotations = set(graph.annotation_pattern) | set(graph.annotation_url) | set(graph.annotation_owl) | set(
        graph.annotation_list)
    return defined_annotations - get_annotations(graph)


def get_unused_list_annotation_values(graph):
    """Gets all of the unused values for list annotations
    
    :param pybel.BELGraph graph: A BEL graph
    :return: A dictionary of {str annotation: set of str values that aren't used}
    :rtype: dict[str, set[str]]
    """
    result = {}
    for annotation, values in graph.annotation_list.items():
        used_values = get_annotation_values(graph, annotation)
        if len(used_values) == len(values):  # all values have been used
            continue
        result[annotation] = set(values) - used_values
    return result


def _get_annotation_values_by_annotation_helper(graph):
    """Gets the list of values, with duplicates, for each annotation used in a BEL graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A dictionary of {annotation key: list of annotation values}
    :rtype: dict[str, list[str]]
    """
    result = defaultdict(set)

    for data in graph_edge_data_iter(graph):
        if ANNOTATIONS not in data:
            continue

        for key, value in data[ANNOTATIONS].items():
            result[key].add(value)

    return dict(result)


def get_annotation_values_by_annotation(graph):
    """Gets the set of values for each annotation used in a BEL graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A dictionary of {annotation key: set of annotation values}
    :rtype: dict[str, set[str]]
    """
    return get_value_sets(_get_annotation_values_by_annotation_helper(graph))


def get_annotations_containing_keyword(graph, keyword):
    """Gets annotation/value pairs for values for whom the search string is a substring

    :param pybel.BELGraph graph: A BEL graph
    :param str keyword: Search for annotations whose values have this as a substring
    :rtype: list[dict[str,str]
    """
    return [
        {
            'annotation': annotation,
            'value': value
        }
        for annotation, values in get_annotation_values_by_annotation(graph).items()
        for value in values
        if keyword.lower() in value.lower()
    ]


def count_annotation_values(graph, annotation):
    """Counts in how many edges each annotation appears in a graph

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to count
    :return: A Counter from {annotation value: frequency}
    :rtype: collections.Counter
    """
    return Counter(
        d[ANNOTATIONS][annotation]
        for _, _, d in graph.edges_iter(data=True)
        if check_has_annotation(d, annotation)
    )


def get_annotation_values(graph, annotation):
    """Get all values for the given annotation

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to summarize
    :return: A set of all annotation values
    :rtype: set[str]
    """
    return set(count_annotation_values(graph, annotation))


def count_annotation_values_filtered(graph, annotation, source_filter=None, target_filter=None):
    """Counts in how many edges each annotation appears in a graph, but filter out source nodes and target nodes

    See :func:`pybel_tools.utils.keep_node` for a basic filter.

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to count
    :param source_filter: A predicate (graph, node) -> bool for keeping source nodes
    :type source_filter: types.FunctionType
    :param target_filter: A predicate (graph, node) -> bool for keeping target nodes
    :type target_filter: types.FunctionType
    :return: A Counter from {annotation value: frequency}
    :rtype: Counter
    """
    source_filter = keep_node_permissive if source_filter is None else source_filter
    target_filter = keep_node_permissive if target_filter is None else target_filter

    return Counter(
        d[ANNOTATIONS][annotation]
        for u, v, d in graph.edges_iter(data=True)
        if check_has_annotation(d, annotation) and source_filter(graph, u) and target_filter(graph, v)
    )


def _iter_pairs(graph):
    """Iterates over unique node-node pairs in the graph
    
    :param pybel.BELGraph graph: A BEL graph
    :rtype: iter
    """
    for u, v in set(graph.edges_iter()):
        yield u, v


def get_all_relations(graph, u, v):
    """Returns the set of all relations between a given pair of nodes
    
    :param pybel.BELGraph graph: A BEL graph
    :param tuple u: The source BEL node
    :param tuple v: The target BEL node
    :rtype: set[str]
    """
    return {
        data[RELATION]
        for data in graph.edge[u][v].values()
    }


def pair_is_consistent(graph, u, v):
    """Returns if the edges between the given nodes are consistent, meaning they all have the same relation

    :param pybel.BELGraph graph: A BEL graph
    :param tuple u: The source BEL node
    :param tuple v: The target BEL node
    :return: If the edges aren't consistent, return false, otherwise return the relation type
    :rtype: bool or str
    """
    relations = get_all_relations(graph, u, v)

    if 1 != len(relations):
        return False

    return list(relations)[0]


def relation_set_has_contradictions(relations):
    """Returns if the set of relations contains a contradiction

    :param set[str] relations: A set of relations
    :rtype: bool
    """
    has_increases = any(relation in CAUSAL_INCREASE_RELATIONS for relation in relations)
    has_decreases = any(relation in CAUSAL_DECREASE_RELATIONS for relation in relations)
    has_cnc = any(relation == CAUSES_NO_CHANGE for relation in relations)
    return 1 < sum([has_cnc, has_decreases, has_increases])


def pair_has_contradiction(graph, u, v):
    """Checks if a pair of nodes has any contradictions in their causal relationships.
    
    :param pybel.BELGraph graph: A BEL graph
    :param tuple u: The source BEL node
    :param tuple v: The target BEL node
    :return: Do the edges between these nodes have a contradiction?
    :rtype: bool
    """
    relations = get_all_relations(graph, u, v)
    return relation_set_has_contradictions(relations)


def get_contradictory_pairs(graph):
    """Iterates over contradictory node pairs in the graph based on their causal relationships
    
    :param pybel.BELGraph graph: A BEL graph
    :return: An iterator over (source, target) node pairs that have contradictory causal edges
    :rtype: iter
    """
    for u, v in _iter_pairs(graph):
        if pair_has_contradiction(graph, u, v):
            yield u, v


def get_consistent_edges(graph):
    """Yields pairs of (source node, target node) for which all of their edges have the same type of relation.

    :param pybel.BELGraph graph: A BEL graph
    :return: An iterator over (source, target) node pairs corresponding to edges with many inconsistent relations
    :rtype: iter[tuple]
    """
    for u, v in _iter_pairs(graph):
        if pair_is_consistent(graph, u, v):
            yield u, v


def get_inconsistent_edges(graph):
    """Returns an iterator over inconsistent edges

    :param pybel.BELGraph graph: A BEL graph
    :return: An iterator over (source, target) node pairs corresponding to edges with many inconsistent relations
    :rtype: iter
    """
    for u, v in _iter_pairs(graph):
        if not pair_is_consistent(graph, u, v):
            yield u, v


def _pathology_iterator(graph):
    """Iterates over the diseases encountered in edges
    
    :param pybel.BELGraph graph: A BEL graph
    :rtype: iter
    """
    for u, v in _iter_pairs(graph):
        if graph.node[u][FUNCTION] == PATHOLOGY:
            yield u
        if graph.node[v][FUNCTION] == PATHOLOGY:
            yield v


def count_pathologies(graph):
    """Returns a counter of all of the mentions of pathologies in a network

    :param pybel.BELGraph graph: A BEL graph
    :rtype: Counter
    """
    return Counter(_pathology_iterator(graph))


def get_tree_annotations(graph):
    """Builds tree structure with annotation for a given graph
    
    :param pybel.BELGraph graph: A BEL Graph
    :return: The JSON structure necessary for building the tree box
    :rtype: list[dict]
    """
    annotations = get_annotation_values_by_annotation(graph)
    return [
        {'text': annotation, 'children': [{'text': value} for value in sorted(values)]}
        for annotation, values in sorted(annotations.items())
    ]
