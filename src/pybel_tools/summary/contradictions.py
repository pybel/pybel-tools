# -*- coding: utf-8 -*-

"""Functions for identifying contradictions."""

from typing import Collection

from pybel import BELGraph
from pybel.constants import CAUSAL_DECREASE_RELATIONS, CAUSAL_INCREASE_RELATIONS, CAUSES_NO_CHANGE, RELATION
from pybel.dsl import BaseEntity

__all__ = [
    'pair_has_contradiction',
    'relation_set_has_contradictions',
]


def pair_has_contradiction(graph: BELGraph, u: BaseEntity, v: BaseEntity) -> bool:
    """Check if a pair of nodes has any contradictions in their causal relationships.

    Assumes both nodes are in the graph.
    """
    relations = {data[RELATION] for data in graph[u][v].values()}
    return relation_set_has_contradictions(relations)


# TODO should this consider correlations?
def relation_set_has_contradictions(relations: Collection[str]) -> bool:
    """Return if the set of BEL relations contains a contradiction."""
    has_increases = any(relation in CAUSAL_INCREASE_RELATIONS for relation in relations)
    has_decreases = any(relation in CAUSAL_DECREASE_RELATIONS for relation in relations)
    has_cnc = any(relation == CAUSES_NO_CHANGE for relation in relations)
    return 1 < sum([has_cnc, has_decreases, has_increases])
