# -*- coding: utf-8 -*-

"""Performs concordance analysis."""

from __future__ import annotations

import enum
import logging
from collections import Counter
from dataclasses import dataclass
from functools import partial
from typing import List, Mapping, Optional, Tuple

from pybel import BELGraph, BaseEntity
from pybel.constants import (
    CAUSAL_DECREASE_RELATIONS, CAUSAL_INCREASE_RELATIONS, CAUSES_NO_CHANGE, NEGATIVE_CORRELATION, POSITIVE_CORRELATION,
    RELATION,
)
from pybel.struct import get_subgraphs_by_annotation
from pybel.struct.mutation import collapse_all_variants, collapse_to_genes
from ..mutation.random import random_by_edges, shuffle_node_data, shuffle_relations

__all__ = [
    'Concordance',
    'edge_concords',
    'calculate_concordance_helper',
    'calculate_concordance',
    'calculate_concordance_by_annotation',
    'calculate_concordance_probability',
    'calculate_concordance_probability_by_annotation',
]

log = logging.getLogger(__name__)

UP = CAUSAL_INCREASE_RELATIONS | {POSITIVE_CORRELATION}
DOWN = CAUSAL_DECREASE_RELATIONS | {NEGATIVE_CORRELATION}


def get_cutoff(value: float, cutoff: Optional[float] = None) -> int:
    """Assign if a value is greater than or less than a cutoff."""
    cutoff = cutoff if cutoff is not None else 0

    if value > cutoff:
        return 1

    if value < (-1 * cutoff):
        return -1

    return 0


class Concordance(enum.Enum):
    """Represents the possible results from differential gene expression edge concordance analysis."""

    correct = 0
    incorrect = 1
    ambiguous = 2
    unassigned = 3


@dataclass(frozen=True)
class ConcordanceResult:
    """Stores the results of the concordance."""

    correct: int
    incorrect: int
    ambiguous: int
    unassigned: int

    @classmethod
    def from_iterable(cls, it) -> ConcordanceResult:
        """Generate a result package from an iterable."""
        return cls.from_counter(Counter(it))

    @staticmethod
    def from_counter(c: Counter) -> ConcordanceResult:
        """Generate a result package from a given counter."""
        return ConcordanceResult(
            c[Concordance.correct],
            c[Concordance.incorrect],
            c[Concordance.ambiguous],
            c[Concordance.unassigned],
        )


def edge_concords(
        graph: BELGraph,
        u: BaseEntity,
        v: BaseEntity,
        k: str,
        key: str,
        cutoff: Optional[float] = None,
) -> Concordance:
    """Calculate the concordance of a given edge with the data in the given key.

    :param graph: A BEL graph
    :param u:
    :param v:
    :param k:
    :param key: The node data dictionary key storing the logFC
    :param cutoff: The optional logFC cutoff for significance
    """
    if key not in graph.nodes[u] or key not in graph.nodes[v]:
        return Concordance.unassigned

    relation = graph[u][v][k][RELATION]

    if relation not in (UP | DOWN | {CAUSES_NO_CHANGE}):
        return Concordance.unassigned

    source_regulation = get_cutoff(graph.nodes[u][key], cutoff=cutoff)
    target_regulation = get_cutoff(graph.nodes[v][key], cutoff=cutoff)

    if source_regulation == 1:
        if target_regulation == 1 and relation in UP:
            return Concordance.correct

        elif target_regulation == 1 and relation in DOWN:
            return Concordance.incorrect

        elif target_regulation == -1 and relation in UP:
            return Concordance.incorrect

        elif target_regulation == -1 and relation in DOWN:
            return Concordance.correct

        elif target_regulation == 0 and relation in (DOWN | UP):
            return Concordance.incorrect

        elif target_regulation == 0 and relation == CAUSES_NO_CHANGE:
            return Concordance.correct

        elif target_regulation in {1, -1} and relation == CAUSES_NO_CHANGE:
            return Concordance.incorrect

        else:
            log.warning('%s %s %s %s %s', u, source_regulation, relation, v, target_regulation)
            return Concordance.ambiguous

    elif source_regulation == -1:
        if target_regulation == 1 and relation in UP:
            return Concordance.incorrect

        elif target_regulation == 1 and relation in DOWN:
            return Concordance.correct

        elif target_regulation == -1 and relation in UP:
            return Concordance.correct

        elif target_regulation == -1 and relation in DOWN:
            return Concordance.incorrect

        elif target_regulation == 0 and relation in (DOWN | UP):
            return Concordance.incorrect

        elif target_regulation == 0 and relation == CAUSES_NO_CHANGE:
            return Concordance.correct

        elif target_regulation in {1, -1} and relation == CAUSES_NO_CHANGE:
            return Concordance.incorrect

        else:
            log.warning('%s %s %s %s %s', u, source_regulation, relation, v, target_regulation)
            return Concordance.ambiguous

    else:  # source_regulation == 0
        if target_regulation == 0 and relation == CAUSES_NO_CHANGE:
            return Concordance.correct

        return Concordance.ambiguous


def calculate_concordance_helper(
        graph: BELGraph,
        key: str,
        cutoff: Optional[float] = None,
) -> ConcordanceResult:
    """Help calculate network-wide concordance.

    Assumes data already annotated with given key

    :param graph: A BEL graph
    :param key: The node data dictionary key storing the logFC
    :param cutoff: The optional logFC cutoff for significance
    """
    return ConcordanceResult.from_iterable(
        edge_concords(graph, u, v, k, key, cutoff=cutoff)
        for u, v, k, d in graph.edges(keys=True, data=True)
    )


def calculate_concordance(
        graph: BELGraph,
        key: str,
        cutoff: Optional[float] = None,
        use_ambiguous: bool = False,
) -> float:
    """Calculate the network-wide concordance.

    Assumes data already annotated with given key

    :param graph: A BEL graph
    :param key: The node data dictionary key storing the logFC
    :param cutoff: The optional logFC cutoff for significance
    :param use_ambiguous: Compare to ambiguous edges as well
    """
    correct, incorrect, ambiguous, _ = calculate_concordance_helper(graph, key, cutoff=cutoff)

    try:
        return correct / (correct + incorrect + (ambiguous if use_ambiguous else 0))
    except ZeroDivisionError:
        return -1.0


def one_sided(value: float, distribution: List[float]) -> float:
    """Calculate the one-sided probability of getting a value more extreme than the distribution."""
    assert distribution
    return sum(value < element for element in distribution) / len(distribution)


ConcordanceTest = Tuple[float, List[float], float]


def calculate_concordance_probability(
        graph: BELGraph,
        key: str,
        cutoff: Optional[float] = None,
        permutations: Optional[int] = None,
        percentage: Optional[float] = None,
        use_ambiguous: bool = False,
        permute_type: Optional[str] = None,
) -> ConcordanceTest:
    """Calculate a graph's concordance as well as its statistical probability.

    :param graph: A BEL graph
    :param key: The node data dictionary key storing the logFC
    :param cutoff: The optional logFC cutoff for significance
    :param permutations: The number of random permutations to test. Defaults to 500
    :param percentage: The percentage of the graph's edges to maintain. Defaults to 0.9
    :param use_ambiguous: Compare to ambiguous edges as well
    :param permute_type: Which permutation algorithm should be used
    :returns: A triple of the concordance score, the null distribution, and the p-value.
    """
    if permute_type == 'random_by_edges':
        permute_func = partial(random_by_edges, percentage=percentage)
    elif permute_type == 'shuffle_node_data' or permute_type is None:
        permute_func = partial(shuffle_node_data, key=key, percentage=percentage)
    elif permute_type == 'shuffle_relations':
        permute_func = partial(shuffle_relations, percentage=percentage)
    else:
        raise ValueError('Invalid permute_type: {}'.format(permute_type))

    graph: BELGraph = graph.copy()
    collapse_to_genes(graph)
    collapse_all_variants(graph)

    score = calculate_concordance(graph, key, cutoff=cutoff)

    null_distribution = [
        calculate_concordance(
            graph=permute_func(graph),
            key=key,
            cutoff=cutoff,
            use_ambiguous=use_ambiguous,
        )
        for _ in range(permutations or 500)
    ]

    one_sided_score = one_sided(score, null_distribution)

    return score, null_distribution, one_sided_score


def calculate_concordance_by_annotation(
        graph: BELGraph,
        annotation: str,
        key: str,
        cutoff: Optional[float] = None,
) -> Mapping[str, float]:
    """Return the concordance scores for each stratified graph based on the given annotation.

    :param graph: A BEL graph
    :param annotation: The annotation to group by.
    :param key: The node data dictionary key storing the logFC
    :param cutoff: The optional logFC cutoff for significance
    """
    x = get_subgraphs_by_annotation(graph, annotation)
    return {
        value: calculate_concordance(
            subgraph,
            key,
            cutoff=cutoff,
        )
        for value, subgraph in x.items()
    }


# TODO multithread this
def calculate_concordance_probability_by_annotation(
        graph: BELGraph,
        annotation: str,
        key: str,
        cutoff: Optional[float] = None,
        permutations: Optional[int] = None,
        percentage: Optional[float] = None,
        use_ambiguous: bool = False,
) -> Mapping[str, ConcordanceTest]:
    """Return the results of concordance analysis on each subgraph, stratified by the given annotation.

    :param graph: A BEL graph
    :param annotation: The annotation to group by.
    :param key: The node data dictionary key storing the logFC
    :param cutoff: The optional logFC cutoff for significance
    :param permutations: The number of random permutations to test. Defaults to 500
    :param percentage: The percentage of the graph's edges to maintain. Defaults to 0.9
    :param use_ambiguous: Compare to ambiguous edges as well
    """
    x = get_subgraphs_by_annotation(graph, annotation)
    return {
        value: calculate_concordance_probability(
            subgraph,
            key,
            cutoff=cutoff,
            permutations=permutations,
            percentage=percentage,
            use_ambiguous=use_ambiguous,
        )
        for value, subgraph in x.items()
    }
