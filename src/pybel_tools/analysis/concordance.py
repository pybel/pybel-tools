# -*- coding: utf-8 -*-

"""Performs concordance analysis"""

import enum
import logging
from collections import defaultdict
from functools import partial

from pybel.constants import (
    RELATION,
    CAUSAL_DECREASE_RELATIONS,
    CAUSAL_INCREASE_RELATIONS,
    CAUSES_NO_CHANGE,
    POSITIVE_CORRELATION,
    NEGATIVE_CORRELATION,
)
from ..mutation import collapse_all_variants, collapse_by_central_dogma_to_genes_out_place
from ..mutation.random import random_by_edges, shuffle_relations, shuffle_node_data
from ..selection import get_subgraphs_by_annotation

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


def get_cutoff(value, cutoff=None):
    """Assigns if a value is greater than or less than a cutoff

    :param float value:
    :param float cutoff:
    :rtype: int
    """
    cutoff = cutoff if cutoff is not None else 0

    if value > cutoff:
        return 1

    if value < (-1 * cutoff):
        return - 1

    return 0


class Concordance(enum.Enum):
    """Represents the possible results from differential gene expression edge concordance analysis"""
    correct = 0
    incorrect = 1
    ambiguous = 2
    unassigned = 3


def edge_concords(graph, u, v, k, d, key, cutoff=None):
    """

    :param pybel.BELGraph graph: A BEL graph
    :param u:
    :param v:
    :param k:
    :param d:
    :param str key: The node data dictionary key storing the logFC
    :param float cutoff: The optional logFC cutoff for significance
    :rtype: Concordance
    """
    if key not in graph.node[u] or key not in graph.node[v]:
        return Concordance.unassigned

    relation = d[RELATION]

    if relation not in (UP | DOWN | {CAUSES_NO_CHANGE}):
        return Concordance.unassigned

    source_regulation = get_cutoff(graph.node[u][key], cutoff=cutoff)
    target_regulation = get_cutoff(graph.node[v][key], cutoff=cutoff)

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


def calculate_concordance_helper(graph, key, cutoff=None):
    """Helps calculate network-wide concordance

    Assumes data already annotated with given key

    :param pybel.BELGraph graph: A BEL graph
    :param str key: The node data dictionary key storing the logFC
    :param float cutoff: The optional logFC cutoff for significance
    :rtype: tuple[int]
    """
    scores = defaultdict(int)

    for u, v, k, d in graph.edges_iter(keys=True, data=True):
        c = edge_concords(graph, u, v, k, d, key, cutoff=cutoff)
        scores[c] += 1

    return (
        scores[Concordance.correct],
        scores[Concordance.incorrect],
        scores[Concordance.ambiguous],
        scores[Concordance.unassigned],
    )


def calculate_concordance(graph, key, cutoff=None, use_ambiguous=False):
    """Calculates network-wide concordance.

    Assumes data already annotated with given key

    :param pybel.BELGraph graph: A BEL graph
    :param str key: The node data dictionary key storing the logFC
    :param float cutoff: The optional logFC cutoff for significance
    :param bool use_ambiguous: Compare to ambiguous edges as well
    :rtype: float
    """
    correct, incorrect, ambiguous, _ = calculate_concordance_helper(graph, key, cutoff=cutoff)

    try:
        return correct / (correct + incorrect + (ambiguous if use_ambiguous else 0))
    except ZeroDivisionError:
        return -1.0


def one_sided(value, distribution):
    """Calculates the one-sided probability of getting a value more extreme than the distribution

    :param float value:
    :param list[float] distribution:
    :rtype: float
    """
    assert distribution
    return sum(value < element for element in distribution) / len(distribution)


def calculate_concordance_probability(graph, key, cutoff=None, permutations=None, percentage=None, use_ambiguous=False,
                                      permute_type='shuffle_node_data'):
    """Calculates a graph's concordance as well as its statistical probability

    :param pybel.BELGraph graph: A BEL graph
    :param str key: The node data dictionary key storing the logFC
    :param float cutoff: The optional logFC cutoff for significance
    :param int permutations: The number of random permutations to test. Defaults to 500
    :param float percentage: The percentage of the graph's edges to maintain. Defaults to 0.9
    :param bool use_ambiguous: Compare to ambiguous edges as well
    :rtype: tuple
    """
    if permute_type == 'random_by_edges':
        permute_func = partial(random_by_edges, percentage=percentage)
    elif permute_type == 'shuffle_node_data':
        permute_func = partial(shuffle_node_data, key=key, percentage=percentage)
    elif permute_type == 'shuffle_relations':
        permute_func = partial(shuffle_relations, percentage=percentage)
    else:
        raise ValueError

    graph = collapse_by_central_dogma_to_genes_out_place(graph)
    collapse_all_variants(graph)

    score = calculate_concordance(graph, key, cutoff=cutoff)

    distribution = []

    for _ in range(permutations or 500):
        permuted_graph = permute_func(graph)
        permuted_graph_scores = calculate_concordance(permuted_graph, key, cutoff=cutoff, use_ambiguous=use_ambiguous)
        distribution.append(permuted_graph_scores)

    return score, distribution, one_sided(score, distribution)


def calculate_concordance_by_annotation(graph, annotation, key, cutoff=None):
    """Returns the concordance scores for each stratified graph based on the given annotation

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to group by.
    :param str key: The node data dictionary key storing the logFC
    :param float cutoff: The optional logFC cutoff for significance
    :rtype: dict[str,tuple]
    """
    return {
        value: calculate_concordance(subgraph, key, cutoff=cutoff)
        for value, subgraph in get_subgraphs_by_annotation(graph, annotation).items()
    }


# TODO multithread this
def calculate_concordance_probability_by_annotation(graph, annotation, key, cutoff=None, permutations=None,
                                                    percentage=None,
                                                    use_ambiguous=False):
    """Returns the results of concordance analysis on each subgraph, stratified by the given annotation.

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to group by.
    :param str key: The node data dictionary key storing the logFC
    :param float cutoff: The optional logFC cutoff for significance
    :param int permutations: The number of random permutations to test. Defaults to 500
    :param float percentage: The percentage of the graph's edges to maintain. Defaults to 0.9
    :param bool use_ambiguous: Compare to ambiguous edges as well
    :rtype: dict[str,tuple]
    """
    result = [
        (value, calculate_concordance_probability(
            subgraph,
            key,
            cutoff=cutoff,
            permutations=permutations,
            percentage=percentage,
            use_ambiguous=use_ambiguous,
        ))
        for value, subgraph in get_subgraphs_by_annotation(graph, annotation).items()
    ]

    return dict(result)
