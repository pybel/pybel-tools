# -*- coding: utf-8 -*-

"""An implementation of the CausalR algorithm described by [Bradley2017]_.

.. [Bradley2017] Bradley, G., & Barrett, S. J. (2017). `CausalR - extracting mechanistic sense from genome scale
                 data <https://doi.org/10.1093/bioinformatics/btx425>`_. Bioinformatics, (June), 1–3.
"""

from __future__ import annotations

import enum
import logging
from functools import reduce
from operator import itemgetter
from typing import Dict, Mapping, Tuple

import networkx as nx

from pybel import BELGraph, BaseEntity
from pybel.constants import (
    ANALOGOUS_TO, ASSOCIATION, BIOMARKER_FOR, CAUSES_NO_CHANGE, DECREASES, DIRECTLY_DECREASES,
    DIRECTLY_INCREASES, EQUIVALENT_TO, HAS_COMPONENT, HAS_MEMBER, HAS_PRODUCT, HAS_REACTANT, HAS_VARIANT, INCREASES,
    IS_A, NEGATIVE_CORRELATION, POSITIVE_CORRELATION, PROGONSTIC_BIOMARKER_FOR, RATE_LIMITING_STEP_OF, REGULATES,
    RELATION, SUBPROCESS_OF, TRANSCRIBED_TO, TRANSLATED_TO,
)
from pybel.typing import EdgeData
from ..summary.contradictions import pair_has_contradiction
from ..utils import pairwise

log = logging.getLogger(__name__)

__all__ = [
    'rank_causalr_hypothesis',
    'run_cna',
    'get_path_effect',
    'rank_edges',
]

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
    CAUSES_NO_CHANGE: 0,
    REGULATES: 0,
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
    EQUIVALENT_TO: 0,
}


class Effect(enum.Enum):
    """Represents the possible effect of a root node on a target."""

    inhibition = -1
    no_effect = 0
    activation = 1
    ambiguous = None

    @classmethod
    def is_causal(cls, effect: Effect) -> bool:
        """Check if the effect is either inhibition or activation."""
        return effect is cls.inhibition or effect is cls.activation


# TODO combine with data integration module so you can give a graph and a path to a data file

HypothesisMapping = Dict[str, int]


def rank_causalr_hypothesis(
    graph: BELGraph,
    node_to_regulation: Mapping[BaseEntity, int],
    regulator_node: BaseEntity,
) -> Tuple[HypothesisMapping, HypothesisMapping]:
    """Test the regulator hypothesis of the given node on the input data using the algorithm.

    Note: this method returns both +/- signed hypotheses evaluated

    Algorithm:

    1. Calculate the shortest path between the regulator node and each node in observed_regulation
    2. Calculate the concordance of the causal network and the observed regulation when there is path
       between target node and regulator node

    :param graph: A causal graph
    :param node_to_regulation: Nodes to score (1,-1,0)
    :param regulator_node:
    :return Dictionaries with hypothesis results (keys: score, correct, incorrect, ambiguous)
    """
    upregulation_hypothesis: HypothesisMapping = {
        'correct': 0,
        'incorrect': 0,
        'ambiguous': 0,
    }
    downregulation_hypothesis: HypothesisMapping = {
        'correct': 0,
        'incorrect': 0,
        'ambiguous': 0,
    }

    targets = [
        node
        for node in node_to_regulation
        if node != regulator_node
    ]

    predicted_regulations = run_cna(graph, regulator_node, targets)  # + signed hypothesis

    for _, target_node, predicted_regulation in predicted_regulations:
        if Effect.is_causal(predicted_regulation) and predicted_regulation.value == node_to_regulation[target_node]:
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


def run_cna(graph: BELGraph, root: BaseEntity, targets, relationship_dict=None):
    """Return the effect from the root to the target nodes represented as {-1, 1}.

    :param graph: A BEL graph
    :param root: The root node
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


def get_path_effect(graph: BELGraph, path, relationship_dict) -> Effect:
    """Calculate the final effect of the root node to the sink node in the path.

    :param graph: A BEL graph
    :param list path: Path from root to sink node
    :param dict relationship_dict: dictionary with relationship effects
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


def rank_edges(edges: Mapping[str, EdgeData], edge_ranking=None):
    """Return the highest ranked edge from a multiedge.

    :param edges: dictionary with all edges between two nodes
    :param dict edge_ranking: A dictionary of {relationship: score}
    :return: Highest ranked edge
    :rtype: tuple: (edge id, relation, score given ranking)
    """
    if edge_ranking is None:
        edge_ranking = default_edge_ranking

    return max(
        (
            (edge_id, edge_data[RELATION], edge_ranking[edge_data[RELATION]])
            for edge_id, edge_data in edges.items()
        ),
        key=itemgetter(2),
    )
