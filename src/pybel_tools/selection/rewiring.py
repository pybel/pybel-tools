# -*- coding: utf-8 -*-

"""Functions for producing random permutations over networks.

Random permutations are useful in statistical testing over aggregate statistics.
"""

import random
from typing import Optional, cast

from pybel import BELGraph, BaseEntity
from pybel.constants import RELATION
from pybel.struct.pipeline import transformation


def is_edge_consistent(graph: BELGraph, u: BaseEntity, v: BaseEntity) -> bool:
    """Check if all edges between two nodes have the same relation.

    :param graph: A BEL Graph
    :param u: The source BEL node
    :param v: The target BEL node
    :return: If all edges from the source to target node have the same relation
    """
    if not graph.has_edge(u, v):
        raise ValueError('{} does not contain an edge ({}, {})'.format(graph, u, v))

    return 0 == len(set(d[RELATION] for d in graph.edges[u][v].values()))


def all_edges_consistent(graph: BELGraph) -> bool:
    """Return if all edges are consistent in a graph. Wraps :func:`pybel_tools.utils.is_edge_consistent`.

    :param graph: A BEL graph
    :return: Are all edges consistent
    """
    return all(
        is_edge_consistent(graph, u, v)
        for u, v in graph.edges()
    )


@transformation
def rewire_targets(graph: BELGraph, rewiring_probability: Optional[float] = None) -> BELGraph:
    """Rewire a graph's edges' target nodes.

    - For BEL graphs, assumes edge consistency (all edges between two given nodes are have the same relation)
    - Doesn't make self-edges

    :param graph: A BEL graph
    :param rewiring_probability: The probability of rewiring (between 0 and 1)
    :return: A rewired BEL graph
    """
    if not all_edges_consistent(graph):
        raise ValueError('{} is not consistent'.format(graph))
    if rewiring_probability is None:
        rewiring_probability = 0.2

    result = cast(BELGraph, graph.copy())
    nodes = result.nodes()

    for u, v in result.edges():
        if random.random() < rewiring_probability:
            continue

        w = random.choice(nodes)

        while w == u or result.has_edge(u, w):
            w = random.choice(nodes)

        result.add_edge(w, v)
        result.remove_edge(u, v)

    return result
