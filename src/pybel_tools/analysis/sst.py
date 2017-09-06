# -*- coding: utf-8 -*-

"""

This module provides an implementation of the sampling of spanning trees (SST) algorithm originally published by
Vasilyev, et al. in https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-7-516

"""

import enum
from functools import reduce
from operator import itemgetter
import networkx as nx
from pybel.constants import *

from ..utils import pairwise
from ..summary.edge_summary import pair_has_contradiction

causal_effect_dict = {
    INCREASES: 1,
    DIRECTLY_INCREASES: 1,
    DECREASES: -1,
    DIRECTLY_DECREASES: -1,
    NEGATIVE_CORRELATION: -1,
    POSITIVE_CORRELATION: 1,
}

default_edge_ranking = {
    INCREASES: 3,
    DIRECTLY_INCREASES: 4,
    DECREASES: 3,
    DIRECTLY_DECREASES: 4,
    RATE_LIMITING_STEP_OF: 0,
    CAUSES_NO_CHANGE: 0, REGULATES: 0,
    NEGATIVE_CORRELATION: 2,
    POSITIVE_CORRELATION: 2,
    ASSOCIATION: 1,
    HAS_MEMBER: 0,
    HAS_PRODUCT: 0,
    HAS_COMPONENT: 0,
    HAS_VARIANT: 0,
    HAS_REACTANT: 0,
    TRANSLATED_TO: 0,
    TRANSCRIBED_TO: 0,
    IS_A: 0,
    SUBPROCESS_OF: 0,
    ANALOGOUS_TO: 0,
    BIOMARKER_FOR: 0,
    PROGONSTIC_BIOMARKER_FOR: 0,
    EQUIVALENT_TO: 0
}


def rank_edges(edges, edge_ranking=None):
    """Returns the highest ranked edge from a multiedge

    :param dict edges: dictionary with all edges between two nodes
    :param dict edge_ranking: A dictionary of {relationship: score}
    :return: Highest ranked edge
    :rtype: tuple: (edge id, relation, score given ranking)
    """

    edge_ranking = default_edge_ranking if edge_ranking is None else edge_ranking

    edges_scores = [
        (edge_id, edge_data[RELATION], edge_ranking[edge_data[RELATION]])
        for edge_id, edge_data in edges.items()
    ]

    return max(edges_scores, key=itemgetter(2))


class Effect(enum.Enum):
    """Represents the possible effect of a root node on a target"""
    inhibition = -1
    no_effect = 0
    activation = 1
    ambiguous = None


def get_path_effect(graph, path, relationship_dict):
    """ Calculates the final effect of the root node to the sink node in the path

    :param pybel.BELGraph graph: A BEL graph
    :param list path: Path from root to sink node
    :param dict relationship_dict: dictionary with relationship effects
    :rtype: Effect
    """

    causal_effect = []

    for predecessor, successor in pairwise(path):

        if pair_has_contradiction(graph, predecessor, successor):
            return Effect.ambiguous

        edges = graph.get_edge_data(predecessor, successor)

        edge_key, edge_relation, _ = rank_edges(edges)

        relation = graph[predecessor][successor][edge_key][RELATION]

        # Returns Effect.no_effect if there is a non causal edge in path
        if relation not in relationship_dict or relationship_dict[relation] == 0:
            return Effect.no_effect

        causal_effect.append(relationship_dict[relation])

    final_effect = reduce(lambda x, y: x * y, causal_effect)

    return Effect.activation if final_effect == 1 else Effect.inhibition


def run_cna(graph, root, targets, relationship_dict=None):
    """ Returns the effect from the root to the target nodes represented as {-1,1}

    :param pybel.BELGraph graph: A BEL graph
    :param tuple root: The root node
    :param iter targets: The targets nodes
    :param dict relationship_dict: dictionary with relationship effects
    :return list[tuple]:
    """

    causal_effects = []

    relationship_dict = causal_effect_dict if relationship_dict is None else relationship_dict

    for target in targets:
        try:
            shortest_paths = nx.all_shortest_paths(graph, source=root, target=target)

            effects_in_path = set()

            for shortest_path in shortest_paths:
                effects_in_path.add(get_path_effect(graph, shortest_path, relationship_dict))

            if len(effects_in_path) == 1:
                causal_effects.append((root, target, next(iter(effects_in_path))))  # Append the only predicted effect

            elif Effect.activation in effects_in_path and Effect.inhibition in effects_in_path:
                causal_effects.append((root, target, Effect.ambiguous))

            elif Effect.activation in effects_in_path and Effect.inhibition not in effects_in_path:
                causal_effects.append((root, target, Effect.activation))

            elif Effect.inhibition in effects_in_path and Effect.activation not in effects_in_path:
                causal_effects.append((root, target, Effect.inhibition))

            else:
                log.warning('Exception in set: {}.'.format(effects_in_path))

        except nx.NetworkXNoPath:
            log.warning('No shortest path between: {} and {}.'.format(root, target))

    return causal_effects


def get_random_walk_spanning_tree(graph):
    """Generates a spanning tree from the directed graph using the random walk appropach proposed independently by
    by Broder (1989) and Aldous (1990). It simply generates random walks until all nodes have been covered.

    Algorithm:

    1. Choose a starting vertex s arbitrarily. Set T_V ← {s} and T_E ← ∅.
    2. Do a simple random walk starting at s. Whenever we cross an edge e = {u, v} with v ∈ V ,
        add v to TV and add e to TE.
    3. Stop the random walk when TV = V . Output T = (T_V , T_E) as our spanning tree

    :param networkx.DiGraph graph: The input graph
    :rtype: networkx.DiGraph

    .. seealso::

        - https://math.dartmouth.edu/~pw/math100w13/kothari.pdf
        - http://keyulux.com/pdf/spanning_tree.pdf

    """
    raise NotImplementedError


def rank_causalr_hypothesis(graph, node_to_regulation, regulator_node):
    """Returns the results of testing the regulator hypothesis of the given node on the input data.
    Note: this method returns both +/- signed hypotheses evaluated
    The algorithm is described in G Bradley et al. (2017)

    Algorithm

    1. Calculate the shortest path between the regulator node and each node in observed_regulation

    2. Calculate the concordance of the causal network and the observed regulation when there is path
    between target node and regulator node

    :param networkx.DiGraph graph: A causal graph
    :param dict node_to_regulation: Nodes to score (1,-1,0)
    :param networkx.DiGraph graph: A causal graph
    :return Dictionaries with hypothesis results (keys: score, correct, incorrect, ambiguous)
    :rtype: dict

    .. seealso::

        - https://doi.org/10.1093/bioinformatics/btx425

    """

    upregulation_hypothesis = {
        'correct': 0,
        'incorrect': 0,
        'ambiguous': 0
    }
    downregulation_hypothesis = {
        'correct': 0,
        'incorrect': 0,
        'ambiguous': 0
    }

    targets = [
        node
        for node in node_to_regulation
        if node != regulator_node
    ]

    predicted_regulations = run_cna(graph, regulator_node, targets)  # + signed hypothesis

    for _, target_node, predicted_regulation in predicted_regulations:

        if (predicted_regulation is Effect.inhibition or predicted_regulation is Effect.activation) and (
                    predicted_regulation.value == node_to_regulation[target_node]):

            upregulation_hypothesis['correct'] += 1
            downregulation_hypothesis['incorrect'] += 1

        elif predicted_regulation is Effect.ambiguous:
            upregulation_hypothesis['ambiguous'] += 1
            downregulation_hypothesis['ambiguous'] += 1

        elif predicted_regulation is Effect.no_effect:
            continue

        else:
            downregulation_hypothesis['correct'] += 1
            upregulation_hypothesis['incorrect'] += 1

    upregulation_hypothesis['score'] = upregulation_hypothesis['correct'] - upregulation_hypothesis['incorrect']
    downregulation_hypothesis['score'] = downregulation_hypothesis['correct'] - downregulation_hypothesis['incorrect']

    return upregulation_hypothesis, downregulation_hypothesis
