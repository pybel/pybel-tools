# -*- coding: utf-8 -*-

import random
from typing import Optional

from pybel import BELGraph
from pybel.struct.pipeline import transformation
from pybel.struct.utils import update_node_helper
from ..utils import safe_add_edge

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
    percentage = percentage or 0.9

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
    percentage = percentage or 0.9
    assert 0 < percentage <= 1

    edges = graph.edges(keys=True)
    n = int(graph.number_of_edges() * percentage)

    subedges = random.sample(edges, n)

    rv = graph.fresh_copy()

    for u, v, k in subedges:
        safe_add_edge(rv, u, v, k, graph[u][v][k])

    update_node_helper(graph, rv)

    return rv


@transformation
def shuffle_node_data(graph: BELGraph, key: str, percentage: Optional[float] = None) -> BELGraph:
    """Shuffle the node's data.

    Useful for permutation testing.

    :param graph: A BEL graph
    :param key: The node data dictionary key
    :param percentage: What percentage of possible swaps to make
    """
    percentage = percentage or 0.3
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
    percentage = percentage or 0.3
    assert 0 < percentage <= 1

    n = graph.number_of_edges()
    swaps = int(percentage * n * (n - 1) / 2)

    result: BELGraph = graph.copy()

    edges = result.edges(keys=True)

    for _ in range(swaps):
        (s1, t1, k1), (s2, t2, k2) = random.sample(edges, 2)
        result[s1][t1][k1], result[s2][t2][k2] = result[s2][t2][k2], result[s1][t1][k1]

    return result
