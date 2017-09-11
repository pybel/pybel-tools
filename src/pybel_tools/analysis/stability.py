# -*- coding: utf-8 -*-

import itertools as itt
import logging

from networkx import DiGraph, Graph

from pybel.constants import *
from ..selection import get_causal_subgraph
from ..summary import get_all_relations, relation_set_has_contradictions

__all__ = [
    'get_contradiction_summary',
    'get_regulatory_pairs',
    'get_chaotic_pairs',
    'get_dampened_pairs',
    'get_correlation_graph',
    'get_correlation_triangles',
    'get_separate_unstable_correlation_triples',
    'get_mutually_unstable_correlation_triples',
    'get_triangles',
    'jens_transformation_alpha',
    'jens_transformation_beta',
    'get_jens_unstable',
    'get_increase_mismatch_triplets',
    'get_decrease_mismatch_triplets',
    'get_chaotic_triplets',
    'get_dampened_triplets',
    'summarize_stability',
]

log = logging.getLogger(__name__)


def get_contradiction_summary(graph):
    """Yield triplets of (source node, target node, set of relations) for (source node, target node) pairs
    that have multiple, contradictory relations.

    :param pybel.BELGraph graph: A BEL graph
    :rtype: iter[tuple]
    """
    for u, v in set(graph.edges_iter()):
        relations = get_all_relations(graph, u, v)
        if relation_set_has_contradictions(relations):
            yield u, v, relations


def get_regulatory_pairs(graph):
    """Finds pairs of nodes that have mutual causal edges that are regulating each other such that ``A -> B`` and
    ``B -| A``.

    :param pybel.BELGraph graph: A BEL graph
    :return: A set of pairs of nodes with mutual causal edges
    :rtype: set
    """
    cg = get_causal_subgraph(graph)

    results = set()

    for u, v, d in cg.edges_iter(data=True):
        if d[RELATION] not in CAUSAL_INCREASE_RELATIONS:
            continue

        if cg.has_edge(v, u) and any(dd[RELATION] in CAUSAL_DECREASE_RELATIONS for dd in cg.edge[v][u].values()):
            results.add((u, v))

    return results


def get_chaotic_pairs(graph):
    """Finds pairs of nodes that have mutual causal edges that are increasing each other such that ``A -> B`` and
    ``B -> A``.

    :param pybel.BELGraph graph: A BEL graph
    :return: A set of pairs of nodes with mutual causal edges
    :rtype: set
    """
    cg = get_causal_subgraph(graph)

    results = set()

    for u, v, d in cg.edges_iter(data=True):
        if d[RELATION] not in CAUSAL_INCREASE_RELATIONS:
            continue

        if cg.has_edge(v, u) and any(dd[RELATION] in CAUSAL_INCREASE_RELATIONS for dd in cg.edge[v][u].values()):
            results.add(tuple(sorted([u, v], key=str)))

    return results


def get_dampened_pairs(graph):
    """Finds pairs of nodes that have mutual causal edges that are decreasing each other such that ``A -| B`` and
    ``B -| A``.

    :param pybel.BELGraph graph: A BEL graph
    :return: A set of pairs of nodes with mutual causal edges
    :rtype: set
    """
    cg = get_causal_subgraph(graph)

    results = set()

    for u, v, d in cg.edges_iter(data=True):
        if d[RELATION] not in CAUSAL_DECREASE_RELATIONS:
            continue

        if cg.has_edge(v, u) and any(dd[RELATION] in CAUSAL_DECREASE_RELATIONS for dd in cg.edge[v][u].values()):
            results.add(tuple(sorted([u, v], key=str)))

    return results


def get_correlation_graph(graph):
    """Extracts a graph of only correlative relationships
    
    :param pybel.BELGraph graph: A BEL Graph
    :rtype: networkx.Graph
    """
    result = Graph()

    for u, v, d in graph.edges_iter(data=True):
        if d[RELATION] not in CORRELATIVE_RELATIONS:
            continue

        if not result.has_edge(u, v):
            result.add_edge(u, v, **{d[RELATION]: True})

        elif d[RELATION] not in result.edge[u][v]:
            log.log(5, 'broken correlation relation for %s, %s', u, v)
            result.edge[u][v][d[RELATION]] = True
            result.edge[v][u][d[RELATION]] = True

    return result


def get_correlation_triangles(graph):
    """Returns a set of all triangles pointed by the given node

    :param networkx.Graph graph: A non-directional graph
    :rtype: set[tuple]
    """
    return {
        tuple(sorted([n, u, v], key=str))
        for n in graph.nodes_iter()
        for u, v in itt.combinations(graph.edge[n], 2)
        if graph.has_edge(u, v)
    }


def get_triangles(graph):
    """Gets a set of triples representing the 3-cycles from a directional graph. Each 3-cycle is returned once,
    with nodes in sorted order.

    :param networkx.DiGraph graph: A directional graph
    :rtype: set[tuple]
    """
    results = set()

    for a, b in graph.edges_iter():
        for c in graph.successors(b):
            if graph.has_edge(c, a):
                canonical_ordering = tuple(sorted([a, b, c], key=str))
                results.add(canonical_ordering)

    return results


def get_separate_unstable_correlation_triples(graph):
    """Yields all triples of nodes A, B, C such that ``A positiveCorrelation B``, ``A positiveCorrelation C``, and
    ``B negativeCorrelation C``

    :param pybel.BELGraph graph: A BEL graph
    :return: An iterator over triples of unstable graphs, where the second two are negative
    :rtype: iter[tuple]
    """
    cg = get_correlation_graph(graph)

    for a, b, c in get_correlation_triangles(cg):
        if POSITIVE_CORRELATION in cg.edge[a][b] and POSITIVE_CORRELATION in cg.edge[b][c] and NEGATIVE_CORRELATION in \
                cg.edge[a][c]:
            yield b, a, c
        if POSITIVE_CORRELATION in cg.edge[a][b] and NEGATIVE_CORRELATION in cg.edge[b][c] and POSITIVE_CORRELATION in \
                cg.edge[a][c]:
            yield a, b, c
        if NEGATIVE_CORRELATION in cg.edge[a][b] and POSITIVE_CORRELATION in cg.edge[b][c] and POSITIVE_CORRELATION in \
                cg.edge[a][c]:
            yield c, a, b


def get_mutually_unstable_correlation_triples(graph):
    """Yields all triples of nodes A, B, C such that ``A negativeCorrelation B``, ``B negativeCorrelation C``, and
    ``C negativeCorrelation A``.

    :param pybel.BELGraph graph: A BEL graph
    :rtype: iter[tuple]
    """
    cg = get_correlation_graph(graph)

    for a, b, c in get_correlation_triangles(cg):
        if all(NEGATIVE_CORRELATION in x for x in (cg.edge[a][b], cg.edge[b][c], cg.edge[a][c])):
            yield a, b, c


def jens_transformation_alpha(graph):
    """Applies Jens' transformation (Type 1) to the graph

    1. Induce a subgraph over causal + correlative edges
    2. Transform edges by the following rules:
        - increases => increases
        - decreases => backwards increases
        - positive correlation => two way increases
        - negative correlation => delete

    The resulting graph can be used to search for 3-cycles, which now symbolize unstable triplets where ``A -> B``,
    ``A -| C`` and ``B positiveCorrelation C``.

    :param pybel.BELGraph graph: A BEL graph
    :rtype: networkx.DiGraph
    """
    result = DiGraph()

    for u, v, k, d in graph.edges_iter(keys=True, data=True):
        relation = d[RELATION]

        if relation == POSITIVE_CORRELATION:
            result.add_edge(u, v)
            result.add_edge(v, u)

        elif relation in CAUSAL_INCREASE_RELATIONS:
            result.add_edge(u, v)

        elif relation in CAUSAL_DECREASE_RELATIONS:
            result.add_edge(v, u)

    for node in result.nodes_iter():
        result.node[node].update(graph.node[node])

    return result


def jens_transformation_beta(graph):
    """Applies Jens' Transformation (Type 2) to the graph

    1. Induce a subgraph over causal and correlative relations
    2. Transform edges with the following rules:
        - increases => backwards decreases
        - decreases => decreases
        - positive correlation => delete
        - negative correlation => two way decreases

    The resulting graph can be used to search for 3-cycles, which now symbolize stable triples where ``A -> B``,
    ``A -| C`` and ``B negativeCorrelation C``.

    :param pybel.BELGraph graph: A BEL graph
    :rtype: networkx.DiGraph
    """
    result = DiGraph()

    for u, v, k, d in graph.edges_iter(keys=True, data=True):
        relation = d[RELATION]

        if relation == NEGATIVE_CORRELATION:
            result.add_edge(u, v)
            result.add_edge(v, u)

        elif relation in CAUSAL_INCREASE_RELATIONS:
            result.add_edge(v, u)

        elif relation in CAUSAL_DECREASE_RELATIONS:
            result.add_edge(u, v)

    for node in result.nodes_iter():
        result.node[node].update(graph.node[node])

    return result


def get_jens_unstable(graph):
    """Yields triples of nodes where ``A -> B``, ``A -| C``, and ``C positiveCorrelation A``. Calculated
    efficiently using the Jens Transformation.

    :param pybel.BELGraph graph: A BEL graph
    :return: An iterable of triplets of nodes
    :rtype: iter[tuple]
    """
    r = jens_transformation_alpha(graph)
    return get_triangles(r)


def _get_mismatch_triplets_helper(graph, relation_set):
    """
    :param pybel.BELGraph graph: A BEL graph
    :return: An iterable of mismatch triples
    :rtype iter[tuple]
    """
    for node in graph.nodes_iter():
        children = {
            v
            for _, v, d in graph.out_edges_iter(node, data=True)
            if d[RELATION] in relation_set
        }

        for a, b in itt.combinations(children, 2):
            if b not in graph.edge[a]:
                continue

            if any(d[RELATION] == NEGATIVE_CORRELATION for d in graph.edge[a][b].values()):
                yield node, a, b


def get_increase_mismatch_triplets(graph):
    """Iterates over triples of nodes where ``A -> B``, ``A -> C``, and ``C negativeCorrelation A``.
    
    :param pybel.BELGraph graph: A BEL graph
    :return: An iterable of triplets of nodes
    :rtype: iter[tuple]
    """
    return _get_mismatch_triplets_helper(graph, CAUSAL_INCREASE_RELATIONS)


def get_decrease_mismatch_triplets(graph):
    """Iterates over triplets of nodes where ``A -| B``, ``A -| C``, and ``C negativeCorrelation A``.

    :param pybel.BELGraph graph: A BEL graph
    :return: An iterable of triplets of nodes
    :rtype: iter[tuple]
    """
    return _get_mismatch_triplets_helper(graph, CAUSAL_DECREASE_RELATIONS)


def _get_disregulated_triplets_helper(graph, relation_set):
    """
    :param pybel.BELGraph graph: A BEL graph
    :param set[str] relation_set: A set of relations to keep 
    :return: 
    """
    result = DiGraph()

    for u, v, d in graph.edges_iter(data=True):
        if d[RELATION] in relation_set:
            result.add_edge(u, v)

    for node in result.nodes_iter():
        result.node[node].update(graph.node[node])

    return get_triangles(result)


def get_chaotic_triplets(graph):
    """Iterates over triples of nodes that mutually increase each other, such as when ``A -> B``, ``B -> C``, and
    ``C -> A``.

    :param pybel.BELGraph graph: A BEL graph
    :return: An iterable of triplets of nodes
    :rtype: iter[tuple]
    """
    return _get_disregulated_triplets_helper(graph, CAUSAL_INCREASE_RELATIONS)


def get_dampened_triplets(graph):
    """Iterates over triples of nodes that mutually decreases each other, such as when ``A -| B``,
    ``B -| C``, and ``C -| A``.

    :param pybel.BELGraph graph: A BEL graph
    :return: An iterable of triplets of nodes
    :rtype: iter[tuple]
    """
    return _get_disregulated_triplets_helper(graph, CAUSAL_DECREASE_RELATIONS)


def summarize_stability(graph):
    """Summarize the stability of the graph

    :param pybel.BELGraph graph: A BEL graph
    :rtype: dict
    """

    regulatory_pairs = get_regulatory_pairs(graph)
    chaotic_pairs = get_chaotic_pairs(graph)
    dampened_pairs = get_dampened_pairs(graph)
    contraditory_pairs = get_contradiction_summary(graph)
    separately_unstable_triples = get_separate_unstable_correlation_triples(graph)
    mutually_unstable_triples = get_mutually_unstable_correlation_triples(graph)
    jens_unstable_triples = get_jens_unstable(graph)
    increase_mismatch_triples = get_increase_mismatch_triplets(graph)
    decrease_mismatch_triples = get_decrease_mismatch_triplets(graph)
    chaotic_triples = get_chaotic_triplets(graph)
    dampened_triples = get_dampened_triplets(graph)

    def count_or_len(it):
        return sum(1 for _ in it)

    return {
        'Regulatory Pairs': count_or_len(regulatory_pairs),
        'Chaotic Pairs': count_or_len(chaotic_pairs),
        'Dampened Pairs': count_or_len(dampened_pairs),
        'Contradictory Pairs': count_or_len(contraditory_pairs),
        'Separately Unstable Triples': count_or_len(separately_unstable_triples),
        'Mutually Unstable Triples': count_or_len(mutually_unstable_triples),
        'Jens Unstable Triples': count_or_len(jens_unstable_triples),
        'Increase Mismatch Triples': count_or_len(increase_mismatch_triples),
        'Decrease Mismatch Triples': count_or_len(decrease_mismatch_triples),
        'Chaotic Triples': count_or_len(chaotic_triples),
        'Dampened Triples': count_or_len(dampened_triples)
    }
