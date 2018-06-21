# -*- coding: utf-8 -*-

import random

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
def random_by_nodes(graph, percentage=None):
    """Gets a random graph by inducing over a percentage of the original nodes

    :param pybel.BELGraph graph: A BEL graph
    :param float percentage: The percentage of edges to keep
    :rtype: pybel.BELGraph
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
def random_by_edges(graph, percentage=None):
    """Gets a random graph by keeping a certain percentage of original edges

    :param pybel.BELGraph graph: A BEL graph
    :param float percentage: The percentage of edges to keep
    :rtype: pybel.BELGraph
    """
    percentage = percentage or 0.9

    assert 0 < percentage <= 1

    edges = graph.edges(keys=True)
    n = int(graph.number_of_edges() * percentage)

    subedges = random.sample(edges, n)

    rv = graph.fresh_copy()

    for u, v, k in subedges:
        safe_add_edge(rv, u, v, k, graph.edge[u][v][k])

    update_node_helper(graph, rv)

    return rv


@transformation
def shuffle_node_data(graph, key, percentage=None):
    """Shuffles the node's data. Useful for permutation testing.

    :param pybel.BELGraph graph: A BEL graph
    :param str key: The node data dictionary key
    :param float percentage: What percentage of possible swaps to make
    :rtype: pybel.BELGraph
    """
    percentage = percentage or 0.3

    assert 0 < percentage <= 1

    n = graph.number_of_nodes()
    swaps = int(percentage * n * (n - 1) / 2)

    result = graph.copy()

    for _ in range(swaps):
        s, t = random.sample(result.node, 2)
        result.node[s][key], result.node[t][key] = result.node[t][key], result.node[s][key]

    return result


@transformation
def shuffle_relations(graph, percentage=None):
    """Shuffles the relations. Useful for permutation testing.

    :param pybel.BELGraph graph: A BEL graph
    :param float percentage: What percentage of possible swaps to make
    :rtype: pybel.BELGraph
    """
    percentage = percentage or 0.3

    assert 0 < percentage <= 1

    n = graph.number_of_edges()
    swaps = int(percentage * n * (n - 1) / 2)

    result = graph.copy()

    edges = result.edges(keys=True)

    for _ in range(swaps):
        (s1, t1, k1), (s2, t2, k2) = random.sample(edges, 2)
        result.edge[s1][t1][k1], result.edge[s2][t2][k2] = result.edge[s2][t2][k2], result.edge[s1][t1][k1]

    return result
