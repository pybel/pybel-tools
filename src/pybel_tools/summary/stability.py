# -*- coding: utf-8 -*-

"""Functions for assessing the stability of a network."""

import itertools as itt
import logging
from typing import Iterable, Mapping, Set, Tuple

from networkx import DiGraph, Graph

from pybel import BELGraph
from pybel.constants import (
    CAUSAL_DECREASE_RELATIONS, CAUSAL_INCREASE_RELATIONS, CORRELATIVE_RELATIONS,
    NEGATIVE_CORRELATION, POSITIVE_CORRELATION, RELATION,
)
from pybel.dsl import BaseEntity
from pybel.struct.utils import update_node_helper
from .contradictions import relation_set_has_contradictions
from ..selection import get_causal_subgraph
from ..typing import NodeTriple, SetOfNodePairs, SetOfNodeTriples

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

logger = logging.getLogger(__name__)


def get_contradiction_summary(graph: BELGraph) -> Set[Tuple[BaseEntity, BaseEntity, Tuple[str]]]:
    """Yield triplets of source, target, relations when there are contradictions."""
    return set(_iterate_contradictions(graph))


def _iterate_contradictions(graph) -> Iterable[Tuple[BaseEntity, BaseEntity, Tuple[str]]]:
    for u, v in set(graph.edges()):
        relations = tuple(sorted({data[RELATION] for data in graph[u][v].values()}))
        if relation_set_has_contradictions(relations):
            yield u, v, relations


def get_regulatory_pairs(graph: BELGraph) -> SetOfNodePairs:
    """Find pairs of nodes such that ``A -> B`` and ``B -| A``.

    :return: A set of pairs of nodes with mutual causal edges
    """
    cg = get_causal_subgraph(graph)

    results = set()

    for u, v, d in cg.edges(data=True):
        if d[RELATION] not in CAUSAL_INCREASE_RELATIONS:
            continue

        if cg.has_edge(v, u) and any(dd[RELATION] in CAUSAL_DECREASE_RELATIONS for dd in cg[v][u].values()):
            results.add((u, v))

    return results


def get_chaotic_pairs(graph: BELGraph) -> SetOfNodePairs:
    """Find pairs of nodes of nodes such that ``A -> B`` and ``B -> A``.

    :return: A set of pairs of nodes with mutual causal edges
    """
    cg = get_causal_subgraph(graph)

    results = set()

    for u, v, d in cg.edges(data=True):
        if d[RELATION] not in CAUSAL_INCREASE_RELATIONS:
            continue

        if cg.has_edge(v, u) and any(dd[RELATION] in CAUSAL_INCREASE_RELATIONS for dd in cg[v][u].values()):
            results.add(tuple(sorted([u, v], key=str)))

    return results


def get_dampened_pairs(graph: BELGraph) -> SetOfNodePairs:
    """Find pairs of nodes such that ``A -| B`` and ``B -| A``.

    :return: A set of pairs of nodes with mutual causal edges
    """
    cg = get_causal_subgraph(graph)

    results = set()

    for u, v, d in cg.edges(data=True):
        if (
            d[RELATION] in CAUSAL_DECREASE_RELATIONS and
            cg.has_edge(v, u) and
            any(dd[RELATION] in CAUSAL_DECREASE_RELATIONS for dd in cg[v][u].values())
        ):
            results.add(tuple(sorted([u, v], key=str)))

    return results


def get_correlation_graph(graph: BELGraph) -> Graph:
    """Extract an undirected graph of only correlative relationships."""
    result = Graph()

    for u, v, d in graph.edges(data=True):
        if d[RELATION] not in CORRELATIVE_RELATIONS:
            continue

        if not result.has_edge(u, v):
            result.add_edge(u, v, **{d[RELATION]: True})

        elif d[RELATION] not in result[u][v]:
            logger.log(5, 'broken correlation relation for %s, %s', u, v)
            result[u][v][d[RELATION]] = True
            result[v][u][d[RELATION]] = True

    return result


def get_correlation_triangles(graph: BELGraph) -> SetOfNodeTriples:
    """Return a set of all triangles pointed by the given node."""
    return {
        tuple(sorted([n, u, v], key=str))
        for n in graph
        for u, v in itt.combinations(graph[n], 2)
        if graph.has_edge(u, v)
    }


def get_triangles(graph: DiGraph) -> SetOfNodeTriples:
    """Get a set of triples representing the 3-cycles from a directional graph.

    Each 3-cycle is returned once, with nodes in sorted order.
    """
    return {
        tuple(sorted([a, b, c], key=str))
        for a, b in graph.edges()
        for c in graph.successors(b)
        if graph.has_edge(c, a)
    }


def get_separate_unstable_correlation_triples(graph: BELGraph) -> SetOfNodeTriples:
    """Yield all triples of nodes A, B, C such that ``A pos B``, ``A pos C``, and ``B neg C``.

    :return: An iterator over triples of unstable graphs, where the second two are negative
    """
    cg = get_correlation_graph(graph)
    return set(_iterate_separate_unstable_correlation_triples(cg))


def _iterate_separate_unstable_correlation_triples(cg) -> Iterable[NodeTriple]:
    for a, b, c in get_correlation_triangles(cg):
        if POSITIVE_CORRELATION in cg[a][b] and POSITIVE_CORRELATION in cg[b][c] and NEGATIVE_CORRELATION in \
                cg[a][c]:
            yield b, a, c
        if POSITIVE_CORRELATION in cg[a][b] and NEGATIVE_CORRELATION in cg[b][c] and POSITIVE_CORRELATION in \
                cg[a][c]:
            yield a, b, c
        if NEGATIVE_CORRELATION in cg[a][b] and POSITIVE_CORRELATION in cg[b][c] and POSITIVE_CORRELATION in \
                cg[a][c]:
            yield c, a, b


def get_mutually_unstable_correlation_triples(graph: BELGraph) -> SetOfNodeTriples:
    """Yield triples of nodes (A, B, C) such that ``A neg B``, ``B neg C``, and ``C neg A``."""
    cg = get_correlation_graph(graph)
    return set(_iterate_mutually_unstable_correlation_triples(cg))


def _iterate_mutually_unstable_correlation_triples(cg) -> Iterable[NodeTriple]:
    for a, b, c in get_correlation_triangles(cg):
        if all(NEGATIVE_CORRELATION in x for x in (cg[a][b], cg[b][c], cg[a][c])):
            yield a, b, c


def jens_transformation_alpha(graph: BELGraph) -> DiGraph:
    """Apply Jens' transformation (Type 1) to the graph.

    1. Induce a sub-graph over causal + correlative edges
    2. Transform edges by the following rules:
        - increases => increases
        - decreases => backwards increases
        - positive correlation => two way increases
        - negative correlation => delete

    The resulting graph can be used to search for 3-cycles, which now symbolize unstable triplets where ``A -> B``,
    ``A -| C`` and ``B positiveCorrelation C``.
    """
    result = DiGraph()

    for u, v, d in graph.edges(data=True):
        relation = d[RELATION]

        if relation == POSITIVE_CORRELATION:
            result.add_edge(u, v)
            result.add_edge(v, u)

        elif relation in CAUSAL_INCREASE_RELATIONS:
            result.add_edge(u, v)

        elif relation in CAUSAL_DECREASE_RELATIONS:
            result.add_edge(v, u)

    return result


def jens_transformation_beta(graph: BELGraph) -> DiGraph:
    """Apply Jens' Transformation (Type 2) to the graph.

    1. Induce a sub-graph over causal and correlative relations
    2. Transform edges with the following rules:
        - increases => backwards decreases
        - decreases => decreases
        - positive correlation => delete
        - negative correlation => two way decreases

    The resulting graph can be used to search for 3-cycles, which now symbolize stable triples where ``A -> B``,
    ``A -| C`` and ``B negativeCorrelation C``.
    """
    result = DiGraph()

    for u, v, d in graph.edges(data=True):
        relation = d[RELATION]

        if relation == NEGATIVE_CORRELATION:
            result.add_edge(u, v)
            result.add_edge(v, u)

        elif relation in CAUSAL_INCREASE_RELATIONS:
            result.add_edge(v, u)

        elif relation in CAUSAL_DECREASE_RELATIONS:
            result.add_edge(u, v)

    return result


def get_jens_unstable(graph: BELGraph) -> SetOfNodeTriples:
    """Yield triples of nodes (A, B, C) where ``A -> B``, ``A -| C``, and ``C positiveCorrelation A``.

    Calculated efficiently using the Jens Transformation.
    """
    r = jens_transformation_alpha(graph)
    return get_triangles(r)


def get_increase_mismatch_triplets(graph: BELGraph) -> SetOfNodeTriples:
    """Yield triples of nodes (A, B, C) where ``A -> B``, ``A -> C``, and ``C negativeCorrelation A``."""
    return set(_iterate_mismatch_triplets(graph, CAUSAL_INCREASE_RELATIONS))


def get_decrease_mismatch_triplets(graph: BELGraph) -> SetOfNodeTriples:
    """Yield triples of nodes (A, B, C) where ``A -| B``, ``A -| C``, and ``C negativeCorrelation A``."""
    return set(_iterate_mismatch_triplets(graph, CAUSAL_DECREASE_RELATIONS))


def _iterate_mismatch_triplets(graph: BELGraph, relation_set: Set[str]) -> Iterable[NodeTriple]:
    for node in graph:
        children = {
            target
            for _, target, data in graph.out_edges(node, data=True)
            if data[RELATION] in relation_set
        }

        for a, b in itt.combinations(children, 2):
            if b not in graph[a]:
                continue

            if any(d[RELATION] == NEGATIVE_CORRELATION for d in graph[a][b].values()):
                yield node, a, b


def get_chaotic_triplets(graph: BELGraph) -> SetOfNodeTriples:
    """Yield triples of nodes (A, B, C) such that ``A -> B``, ``B -> C``, and ``C -> A``."""
    return set(_iterate_disregulated_triplets(graph, CAUSAL_INCREASE_RELATIONS))


def get_dampened_triplets(graph: BELGraph) -> SetOfNodeTriples:
    """Yield triples of nodes (A, B, C) such that ``A -| B``, ``B -| C``, and ``C -| A``."""
    return set(_iterate_disregulated_triplets(graph, CAUSAL_DECREASE_RELATIONS))


def _iterate_disregulated_triplets(graph: BELGraph, relation_set: Set[str]) -> Iterable[NodeTriple]:
    result = DiGraph()

    for u, v, d in graph.edges(data=True):
        if d[RELATION] in relation_set:
            result.add_edge(u, v)

    update_node_helper(graph, result)

    for a, b, c in get_triangles(result):
        if a == b == c:
            continue
        yield a, b, c


def summarize_stability(graph: BELGraph) -> Mapping[str, int]:
    """Summarize the stability of the graph."""
    regulatory_pairs = get_regulatory_pairs(graph)
    chaotic_pairs = get_chaotic_pairs(graph)
    dampened_pairs = get_dampened_pairs(graph)
    contradictory_pairs = get_contradiction_summary(graph)
    separately_unstable_triples = get_separate_unstable_correlation_triples(graph)
    mutually_unstable_triples = get_mutually_unstable_correlation_triples(graph)
    jens_unstable_triples = get_jens_unstable(graph)
    increase_mismatch_triples = get_increase_mismatch_triplets(graph)
    decrease_mismatch_triples = get_decrease_mismatch_triplets(graph)
    chaotic_triples = get_chaotic_triplets(graph)
    dampened_triples = get_dampened_triplets(graph)

    return {
        'Regulatory Pairs': _count_or_len(regulatory_pairs),
        'Chaotic Pairs': _count_or_len(chaotic_pairs),
        'Dampened Pairs': _count_or_len(dampened_pairs),
        'Contradictory Pairs': _count_or_len(contradictory_pairs),
        'Separately Unstable Triples': _count_or_len(separately_unstable_triples),
        'Mutually Unstable Triples': _count_or_len(mutually_unstable_triples),
        'Jens Unstable Triples': _count_or_len(jens_unstable_triples),
        'Increase Mismatch Triples': _count_or_len(increase_mismatch_triples),
        'Decrease Mismatch Triples': _count_or_len(decrease_mismatch_triples),
        'Chaotic Triples': _count_or_len(chaotic_triples),
        'Dampened Triples': _count_or_len(dampened_triples),
    }


def _count_or_len(it):
    return sum(1 for _ in it)
