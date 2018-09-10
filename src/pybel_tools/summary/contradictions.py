# -*- coding: utf-8 -*-

from typing import Set

from pybel import BELGraph
from pybel.constants import CAUSAL_DECREASE_RELATIONS, CAUSAL_INCREASE_RELATIONS, CAUSES_NO_CHANGE, RELATION
from pybel.dsl import BaseEntity


def relation_set_has_contradictions(relations: Set[str]) -> bool:
    """Return if the set of relations contains a contradiction.

    :param relations: A set of BEL relations
    """
    has_increases = any(relation in CAUSAL_INCREASE_RELATIONS for relation in relations)
    has_decreases = any(relation in CAUSAL_DECREASE_RELATIONS for relation in relations)
    has_cnc = any(relation == CAUSES_NO_CHANGE for relation in relations)
    return 1 < sum([has_cnc, has_decreases, has_increases])


def pair_has_contradiction(graph: BELGraph, u: BaseEntity, v) -> bool:
    """Check if a pair of nodes has any contradictions in their causal relationships.

    :param graph: A BEL graph
    :param u: The source BEL node
    :param v: The target BEL node
    :return: Do the edges between these nodes have a contradiction?
    """
    relations = {data[RELATION] for data in graph[u][v].values()}
    return relation_set_has_contradictions(relations)
