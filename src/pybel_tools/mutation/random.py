# -*- coding: utf-8 -*-

"""Random graph permutation functions."""

import random
from typing import Optional

from pybel import BELGraph
from pybel.struct.pipeline import transformation
from pybel.struct.utils import update_node_helper

__all__ = [
    'random_by_nodes',
    'random_by_edges',
    'shuffle_node_data',
    'shuffle_relations',
]


@transformation
def random_by_nodes(graph: BELGraph, percentage: Optional[float] = None) -> BELGraph:
    """Get a random graph by inducing over a percentage of the original nodes.

    :param graph: A BEL graph
    :param percentage: The percentage of edges to keep
    """
    if percentage is None:
        percentage = 0.9
    assert 0 < percentage <= 1

    nodes = graph.nodes()
    n = int(len(nodes) * percentage)

    subnodes = random.sample(nodes, n)

    result = graph.subgraph(subnodes)

    update_node_helper(graph, result)

    return result


@transformation
def random_by_edges(graph: BELGraph, percentage: Optional[float] = None) -> BELGraph:
    """Get a random graph by keeping a certain percentage of original edges.

    :param graph: A BEL graph
    :param percentage: What percentage of eges to take
    """
    if percentage is None:
        percentage = 0.9
    assert 0 < percentage <= 1

    number_edges = int(graph.number_of_edges() * percentage)
    rv = graph.fresh_copy()
    rv.add_edges_from(random.sample(graph.edges(keys=True, data=True), number_edges))
    update_node_helper(graph, rv)
    return rv


@transformation
def shuffle_node_data(graph: BELGraph, key: str, percentage: Optional[float] = None) -> BELGraph:
    """Shuffle the graphs' nodes' data.

    Useful for permutation testing. For example, shuffling differential gene expression values.

    :param graph: A BEL graph
    :param key: The node data dictionary key
    :param percentage: What percentage of possible swaps to make
    """
    if percentage is None:
        percentage = 0.3
    assert 0 < percentage <= 1

    n = graph.number_of_nodes()
    swaps = int(percentage * n * (n - 1) / 2)

    result: BELGraph = graph.copy()

    for _ in range(swaps):
        s, t = random.sample(result.node, 2)
        result.nodes[s][key], result.nodes[t][key] = result.nodes[t][key], result.nodes[s][key]

    return result


@transformation
def shuffle_relations(graph: BELGraph, percentage: Optional[str] = None) -> BELGraph:
    """Shuffle the relations.

    Useful for permutation testing.

    :param graph: A BEL graph
    :param percentage: What percentage of possible swaps to make
    """
    if percentage is None:
        percentage = 0.3
    assert 0 < percentage <= 1

    n = graph.number_of_edges()
    swaps = int(percentage * n * (n - 1) / 2)

    rv = graph.copy()

    edges = rv.edges(keys=True)

    for _ in range(swaps):
        (s1, t1, k1), (s2, t2, k2) = random.sample(edges, 2)
        rv[s1][t1][k1], rv[s2][t2][k2] = rv[s2][t2][k2], rv[s1][t1][k1]

    return rv
