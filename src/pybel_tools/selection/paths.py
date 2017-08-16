# -*- coding: utf-8 -*-


import itertools as itt
from itertools import product, tee

import networkx as nx
from networkx import all_shortest_paths

from pybel.constants import *

__all__ = [
    'get_nodes_in_all_shortest_paths',
    'get_shortest_directed_path_between_subgraphs',
    'get_shortest_undirected_path_between_subgraphs',
]

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


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def rank_path(graph, path, edge_ranking=None):
    """Takes in a path (a list of nodes in the graph) and calculates a score

    :param pybel.BELGraph graph: A BEL graph
    :param list[tuple] path: A list of nodes in the path (includes terminal nodes)
    :param dict edge_ranking: A dictionary of {relationship: score}
    :return: The score for the edge
    :rtype: int
    """
    edge_ranking = default_edge_ranking if edge_ranking is None else edge_ranking

    return sum(max(edge_ranking[d[RELATION]] for d in graph.edge[u][v].values()) for u, v in pairwise(path))


def get_nodes_in_all_shortest_paths(graph, nodes, weight=None):
    """Gets all shortest paths from all nodes to all other nodes in the given list and returns the set of all nodes 
    contained in those paths using :func:`networkx.all_shortest_paths`.

    :param pybel.BELGraph graph: A BEL graph
    :param iter[tuple] nodes: The list of nodes to use to use to find all shortest paths
    :param int cutoff:  Depth to stop the search. Only paths of length <= cutoff are returned.
    :param str weight: Edge data key corresponding to the edge weight. If none, uses unweighted search.
    :return: A set of nodes appearing in the shortest paths between nodes in the BEL graph
    :rtype: set

    .. note:: This can be trivially parallelized using :func:`networkx.single_source_shortest_path`
    """
    return {node
            for u, v in product(nodes, repeat=2)
            for path in all_shortest_paths(graph, u, v, weight=weight)
            for node in path}


# TODO consider all shortest paths?
def _get_shortest__path_between_subgraphs_helper(graph, a, b):
    """Calculate the shortest path that occurs between two disconnected subgraphs A and B going through nodes in
    the source graph

    :param graph: A graph
    :type graph: nx.MultiGraph
    :param a: A subgraph of :code:`graph`, disjoint from :code:`b`
    :type a: nx.MultiGraph
    :param b: A subgraph of :code:`graph`, disjoint from :code:`a`
    :type b: nx.MultiGraph
    :return: A list of the shortest paths between the two subgraphs
    :rtype: list
    """
    shortest_paths = []

    for na, nb in itt.product(a.nodes_iter(), b.nodes_iter()):
        a_b_shortest_path = nx.shortest_path(graph, na, nb)
        shortest_paths.append(a_b_shortest_path)

        b_a_shortest_path = nx.shortest_path(graph, nb, na)
        shortest_paths.append(b_a_shortest_path)

    min_len = min(map(len, shortest_paths))
    return [p for p in shortest_paths if len(p) == min_len]


def get_shortest_directed_path_between_subgraphs(graph, a, b):
    """Calculate the shortest path that occurs between two disconnected subgraphs A and B going through nodes in
    the source graph

    :param pybel.BELGraph graph: A BEL graph
    :param a: A subgraph of :code:`graph`, disjoint from :code:`b`
    :type a: pybel.BELGraph
    :param b: A subgraph of :code:`graph`, disjoint from :code:`a`
    :type b: pybel.BELGraph
    :return: A list of the shortest paths between the two subgraphs
    :rtype: list
    """
    return _get_shortest__path_between_subgraphs_helper(graph, a, b)


def get_shortest_undirected_path_between_subgraphs(graph, a, b):
    """Get the shortest path between two disconnected subgraphs A and B, disregarding directionality of edges in graph

    :param pybel.BELGraph graph: A BEL graph
    :param a: A subgraph of :code:`graph`, disjoint from :code:`b`
    :type a: pybel.BELGraph
    :param b: A subgraph of :code:`graph`, disjoint from :code:`a`
    :type b: pybel.BELGraph
    :return: A list of the shortest paths between the two subgraphs
    :rtype: list
    """
    ug = graph.to_undirected()
    return _get_shortest__path_between_subgraphs_helper(ug, a, b)


def find_root_in_path(graph, path_nodes):
    """Find the 'root' of the path -> The node with the lowest out degree, if multiple:
         root is the one with the highest out degree among those with lowest out degree
    
    :param pybel.BELGraph graph: A BEL Graph
    :param list[tuple] path: A list of nodes in their order in a path
    :rtype: pybel.BELGraph 
    :return: graph: graph of the path
    :rtype: tuple
    :return: root node
      """
    path_graph = graph.subgraph(path_nodes)

    # node_in_degree_tuple: list of tuples with (node,in_degree_of_node) in ascending order
    node_in_degree_tuple = sorted([(n, d) for n, d in path_graph.in_degree().items()], key=lambda x: x[1])
    # node_out_degree_tuple: ordered list of tuples with (node,in_degree_of_node) in descending order
    node_out_degree_tuple = sorted([(n, d) for n, d in path_graph.out_degree().items()], key=lambda x: x[1],
                                   reverse=True)

    # In case all have the same in degree it needs to be reference before
    tied_root_index = 0

    # Get index where the min in_degree stops (in case they are duplicates)
    for i in range(0, (len(node_in_degree_tuple) - 1)):
        if node_in_degree_tuple[i][1] < node_in_degree_tuple[i + 1][1]:
            tied_root_index = i
            break

    # If there are multiple nodes with minimum in_degree take the one with max out degree
    # (in case multiple have the same out degree pick one random)
    if tied_root_index != 0:
        root_tuple = max(node_out_degree_tuple[:tied_root_index], key=lambda x: x[1])
    else:
        root_tuple = node_in_degree_tuple[0]

    return path_graph, root_tuple[0]
