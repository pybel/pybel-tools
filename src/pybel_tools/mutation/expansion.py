# -*- coding: utf-8 -*-

"""Deletion functions to supplement :mod:`pybel.struct.mutation.expansion`."""

import logging
from collections import Counter, defaultdict
from typing import Iterable, Union

import itertools as itt

from pybel import BELGraph
from pybel.constants import *
from pybel.dsl import BaseEntity
from pybel.struct.filters import and_edge_predicates, concatenate_node_predicates, get_nodes_by_function
from pybel.struct.filters.edge_predicates import edge_has_annotation, is_causal_relation
from pybel.struct.filters.node_predicates import keep_node_permissive
from pybel.struct.filters.typing import EdgeIterator, EdgePredicates, NodePredicates
from pybel.struct.pipeline import uni_in_place_transformation
from ..utils import safe_add_edge

__all__ = [
    'get_peripheral_successor_edges',
    'get_peripheral_predecessor_edges',
    'count_sources',
    'count_targets',
    'count_possible_successors',
    'count_possible_predecessors',
    'get_subgraph_edges',
    'get_subgraph_peripheral_nodes',
    'expand_periphery',
    'enrich_complexes',
    'enrich_composites',
    'enrich_reactions',
    'enrich_variants',
    'enrich_unqualified',
    'expand_internal',
    'expand_internal_causal',
]

log = logging.getLogger(__name__)


def get_peripheral_successor_edges(graph: BELGraph, subgraph: BELGraph) -> EdgeIterator:
    """Get the set of possible successor edges peripheral to the sub-graph.

    The source nodes in this iterable are all inside the sub-graph, while the targets are outside.
    """
    for u in subgraph:
        for _, v, k in graph.out_edges(u, keys=True):
            if v not in subgraph:
                yield u, v, k


def get_peripheral_predecessor_edges(graph: BELGraph, subgraph: BELGraph) -> EdgeIterator:
    """Get the set of possible predecessor edges peripheral to the sub-graph.

    The target nodes in this iterable are all inside the sub-graph, while the sources are outside.
    """
    for v in subgraph:
        for u, _, k in graph.in_edges(v, keys=True):
            if u not in subgraph:
                yield u, v, k


def count_sources(edge_iter: EdgeIterator) -> Counter:
    """Count the source nodes in an edge iterator with keys and data.

    :return: A counter of source nodes in the iterable
    """
    return Counter(u for u, _, _ in edge_iter)


def count_targets(edge_iter: EdgeIterator) -> Counter:
    """Count the target nodes in an edge iterator with keys and data.

    :return: A counter of target nodes in the iterable
    """
    return Counter(v for _, v, _ in edge_iter)


def count_possible_successors(graph: BELGraph, subgraph: BELGraph) -> Counter:
    """

    :param graph: A BEL graph
    :param subgraph: An iterator of BEL nodes
    :return: A counter of possible successor nodes
    """
    return count_targets(get_peripheral_successor_edges(graph, subgraph))


def count_possible_predecessors(graph: BELGraph, subgraph: BELGraph) -> Counter:
    """

    :param graph: A BEL graph
    :param subgraph: An iterator of BEL nodes
    :return: A counter of possible predecessor nodes
    """
    return count_sources(get_peripheral_predecessor_edges(graph, subgraph))


def get_subgraph_edges(graph: BELGraph,
                       annotation: str,
                       value: str,
                       source_filter=None,
                       target_filter=None,
                       ):
    """Gets all edges from a given subgraph whose source and target nodes pass all of the given filters

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation:  The annotation to search
    :param str value: The annotation value to search by
    :param source_filter: Optional filter for source nodes (graph, node) -> bool
    :param target_filter: Optional filter for target nodes (graph, node) -> bool
    :return: An iterable of (source node, target node, key, data) for all edges that match the annotation/value and
             node filters
    :rtype: iter[tuple]
    """
    if source_filter is None:
        source_filter = keep_node_permissive

    if target_filter is None:
        target_filter = keep_node_permissive

    for u, v, k, data in graph.edges(keys=True, data=True):
        if not edge_has_annotation(data, annotation):
            continue
        if data[ANNOTATIONS][annotation] == value and source_filter(graph, u) and target_filter(graph, v):
            yield u, v, k, data


def get_subgraph_peripheral_nodes(graph: BELGraph,
                                  subgraph: Iterable[BaseEntity],
                                  node_predicates: NodePredicates = None,
                                  edge_predicates: EdgePredicates = None,
                                  ):
    """Get a summary dictionary of all peripheral nodes to a given sub-graph.

    :return: A dictionary of {external node: {'successor': {internal node: list of (key, dict)},
                                            'predecessor': {internal node: list of (key, dict)}}}
    :rtype: dict

    For example, it might be useful to quantify the number of predecessors and successors:

    >>> from pybel.struct.filters import exclude_pathology_filter
    >>> value = 'Blood vessel dilation subgraph'
    >>> sg = get_subgraph_by_annotation_value(graph, annotation='Subgraph', value=value)
    >>> p = get_subgraph_peripheral_nodes(graph, sg, node_predicates=exclude_pathology_filter)
    >>> for node in sorted(p, key=lambda n: len(set(p[n]['successor']) | set(p[n]['predecessor'])), reverse=True):
    >>>     if 1 == len(p[value][node]['successor']) or 1 == len(p[value][node]['predecessor']):
    >>>         continue
    >>>     print(node,
    >>>           len(p[node]['successor']),
    >>>           len(p[node]['predecessor']),
    >>>           len(set(p[node]['successor']) | set(p[node]['predecessor'])))
    """
    node_filter = concatenate_node_predicates(node_predicates=node_predicates)
    edge_filter = and_edge_predicates(edge_predicates=edge_predicates)

    result = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for u, v, k, d in get_peripheral_successor_edges(graph, subgraph):
        if not node_filter(graph, v) or not node_filter(graph, u) or not edge_filter(graph, u, v, k):
            continue
        result[v]['predecessor'][u].append((k, d))

    for u, v, k, d in get_peripheral_predecessor_edges(graph, subgraph):
        if not node_filter(graph, v) or not node_filter(graph, u) or not edge_filter(graph, u, v, k):
            continue
        result[u]['successor'][v].append((k, d))

    return result


@uni_in_place_transformation
def expand_periphery(universe: BELGraph,
                     graph: BELGraph,
                     node_predicates: NodePredicates = None,
                     edge_predicates: EdgePredicates = None,
                     threshold: int = 2,
                     ) -> None:
    """Iterates over all possible edges, peripheral to a given subgraph, that could be added from the given graph.
    Edges could be added if they go to nodes that are involved in relationships that occur with more than the
    threshold (default 2) number of nodes in the subgraph.

    :param universe: The universe of BEL knowledge
    :param graph: The (sub)graph to expand
    :param threshold: Minimum frequency of betweenness occurrence to add a gap node

    A reasonable edge filter to use is :func:`pybel_tools.filters.keep_causal_edges` because this function can allow
    for huge expansions if there happen to be hub nodes.
    """
    nd = get_subgraph_peripheral_nodes(universe, graph, node_predicates=node_predicates,
                                       edge_predicates=edge_predicates)

    for node, dd in nd.items():
        pred_d = dd['predecessor']
        succ_d = dd['successor']

        in_subgraph_connections = set(pred_d) | set(succ_d)

        if threshold > len(in_subgraph_connections):
            continue

        graph.add_node(node, attr_dict=universe[node])

        for u, edges in pred_d.items():
            for k, d in edges:
                safe_add_edge(graph, u, node, k, d)

        for v, edges in succ_d.items():
            for k, d in edges:
                safe_add_edge(graph, node, v, k, d)


@uni_in_place_transformation
def enrich_complexes(graph: BELGraph) -> None:
    """Add all of the members of the complex abundances to the graph."""
    nodes = list(get_nodes_by_function(graph, COMPLEX))
    for u in nodes:
        for v in u.members:
            graph.add_has_component(u, v)


@uni_in_place_transformation
def enrich_composites(graph: BELGraph):
    """Adds all of the members of the composite abundances to the graph."""
    nodes = list(get_nodes_by_function(graph, COMPOSITE))
    for u in nodes:
        for v in u.members:
            graph.add_has_component(u, v)


@uni_in_place_transformation
def enrich_reactions(graph: BELGraph):
    """Adds all of the reactants and products of reactions to the graph."""
    nodes = list(get_nodes_by_function(graph, REACTION))
    for u in nodes:
        for v in u.reactants:
            graph.add_has_reactant(u, v)

        for v in u.products:
            graph.add_has_product(u, v)


@uni_in_place_transformation
def enrich_variants(graph: BELGraph, func: Union[None, str, Iterable[str]] = None):
    """Add the reference nodes for all variants of the given function.

    :param graph: The target BEL graph to enrich
    :param func: The function by which the subject of each triple is filtered. Defaults to the set of protein, rna,
     mirna, and gene.
    """
    if func is None:
        func = {PROTEIN, RNA, MIRNA, GENE}

    nodes = list(get_nodes_by_function(graph, func))
    for u in nodes:
        parent = u.get_parent()

        if parent is None:
            continue

        if parent not in graph:
            graph.add_has_variant(parent, u)


@uni_in_place_transformation
def enrich_unqualified(graph: BELGraph):
    """Enrich the sub-graph with the unqualified edges from the graph.

    The reason you might want to do this is you induce a sub-graph from the original graph based on an annotation
    filter, but the unqualified edges that don't have annotations that most likely connect elements within your graph
    are not included.

    .. seealso::

        This function thinly wraps the successive application of the following functions:

        - :func:`enrich_complexes`
        - :func:`enrich_composites`
        - :func:`enrich_reactions`
        - :func:`enrich_variants`

    Equivalent to:

    >>> enrich_complexes(graph)
    >>> enrich_composites(graph)
    >>> enrich_reactions(graph)
    >>> enrich_variants(graph)
    """
    enrich_complexes(graph)
    enrich_composites(graph)
    enrich_reactions(graph)
    enrich_variants(graph)


# TODO should this bother checking multiple relationship types?
@uni_in_place_transformation
def expand_internal(universe: BELGraph, graph: BELGraph, edge_predicates: EdgePredicates = None) -> None:
    """Edges between entities in the sub-graph that pass the given filters.

    :param universe: The full graph
    :param graph: A sub-graph to find the upstream information
    :param edge_predicates: Optional list of edge filter functions (graph, node, node, key, data) -> bool
    """
    edge_filter = and_edge_predicates(edge_predicates)

    for u, v in itt.product(graph, repeat=2):
        if graph.has_edge(u, v) or not universe.has_edge(u, v):
            continue

        rs = defaultdict(list)
        for key, data in universe[u][v].items():
            if not edge_filter(universe, u, v, key):
                continue

            rs[data[RELATION]].append((key, data))

        if 1 == len(rs):
            relation = list(rs)[0]
            for key, data in rs[relation]:
                graph.add_edge(u, v, key=key, **data)
        else:
            log.debug('Multiple relationship types found between %s and %s', u, v)


@uni_in_place_transformation
def expand_internal_causal(universe: BELGraph, graph: BELGraph) -> None:
    """Add causal edges between entities in the sub-graph.

    Is an extremely thin wrapper around :func:`expand_internal`.

    :param universe: A BEL graph representing the universe of all knowledge
    :param graph: The target BEL graph to enrich with causal relations between contained nodes

    Equivalent to:

    >>> from pybel_tools.mutation import expand_internal
    >>> from pybel.struct.filters.edge_predicates import is_causal_relation
    >>> expand_internal(universe, graph, edge_predicates=is_causal_relation)
    """
    expand_internal(universe, graph, edge_predicates=is_causal_relation)
