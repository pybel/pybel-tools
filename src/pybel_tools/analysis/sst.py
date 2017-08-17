# -*- coding: utf-8 -*-

"""

This module provides an implementation of the sampling of spanning trees (SST) algorithm originally published by
Vasilyev, et al. in https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-7-516

"""

import networkx as nx
import enum
from functools import reduce
from operator import itemgetter

from pybel.constants import *
from ..utils import pairwise

causal_effect_dict = {
    INCREASES: 1,
    DIRECTLY_INCREASES: 1,
    DECREASES: -1,
    DIRECTLY_DECREASES: -1,
    NEGATIVE_CORRELATION: -1,
    POSITIVE_CORRELATION: 1,
}

default_edge_ranking = {
    INCREASES: 2,
    DIRECTLY_INCREASES: 3,
    DECREASES: 2,
    DIRECTLY_DECREASES: 3,
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


def get_path_effect(graph, path, relationship_dict):
    """ Calculates the final effect of the root node to the sink node in the path

    :param pybel.BELGraph graph: A BEL graph
    :param list path: Path from root to sink node
    :param dict relationship_dict: dictionary with relationship effects
    :rtype: Effect
    """

    causal_effect = []

    for predecessor, successor in pairwise(path):

        edges = graph.get_edge_data(predecessor, successor)

        edge_key, edge_relation, _ = rank_edges(edges)

        relation = graph[predecessor][successor][edge_key][RELATION]

        if relation not in relationship_dict or relationship_dict[relation] == 0:
            return Effect.no_effect

        causal_effect.append(relationship_dict[relation])

    return reduce(lambda x, y: x * y, causal_effect)


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
            path = nx.shortest_path(graph, source=root, target=target)

            causal_effects.append((root, target, get_path_effect(graph, path, relationship_dict)))

        except nx.NetworkXNoPath:
            log.warning('No shortest path between: {} and {}.'.format(root, target))

    return causal_effects
