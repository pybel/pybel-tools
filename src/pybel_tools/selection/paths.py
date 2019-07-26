# -*- coding: utf-8 -*-

"""Path selection tools."""

import itertools as itt
from operator import itemgetter
from typing import Any, List, Mapping, Optional, Tuple

import networkx as nx

from pybel import BELGraph, BaseEntity
from pybel.constants import (
    ANALOGOUS_TO, ASSOCIATION, BIOMARKER_FOR, CAUSES_NO_CHANGE, DECREASES, DIRECTLY_DECREASES, DIRECTLY_INCREASES,
    EQUIVALENT_TO, HAS_COMPONENT, HAS_MEMBER, HAS_PRODUCT, HAS_REACTANT, HAS_VARIANT, INCREASES, IS_A,
    NEGATIVE_CORRELATION, POSITIVE_CORRELATION, PROGONSTIC_BIOMARKER_FOR, RATE_LIMITING_STEP_OF, REGULATES, RELATION,
    SUBPROCESS_OF, TRANSCRIBED_TO, TRANSLATED_TO,
)
from pybel.struct.mutation import get_nodes_in_all_shortest_paths
from ..utils import pairwise

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


def rank_path(graph: BELGraph, path: List[BaseEntity], edge_ranking: Optional[Mapping[str, int]] = None) -> int:
    """Score the given path.

    :param graph: A BEL graph
    :param path: A list of nodes in the path (includes terminal nodes)
    :param edge_ranking: A dictionary of {relationship: score}
    :return: The score for the edge
    """
    if edge_ranking is None:
        edge_ranking = default_edge_ranking

    return sum(
        max(
            edge_ranking[data[RELATION]]
            for data in graph.edges[source][target].values()
        )
        for source, target in pairwise(path)
    )


# TODO consider all shortest paths?
def _get_shortest_path_between_subgraphs_helper(
        graph: nx.Graph,
        a: nx.Graph,
        b: nx.Graph,
) -> List[List[Any]]:
    """Calculate the shortest path(s) between disconnected sub-graphs ``a`` and ``b`` through ``graph``.

    :param graph: A graph
    :param a: A sub-graph of :code:`graph`, disjoint from :code:`b`
    :param b: A sub-graph of :code:`graph`, disjoint from :code:`a`
    :return: A list of the shortest paths between the two sub-graphs
    """
    if graph.is_directed():
        shortest_paths = [
            shortest_path
            for na, nb in itt.product(a, b)
            for shortest_path in (  # do it going both ways because it's directed
                nx.shortest_path(graph, na, nb),
                nx.shortest_path(graph, nb, na),
            )
        ]
    else:
        shortest_paths = [
            nx.shortest_path(graph, na, nb)
            for na, nb in itt.product(a, b)
        ]

    min_len = min(map(len, shortest_paths))
    return [
        shortest_path
        for shortest_path in shortest_paths
        if len(shortest_path) == min_len
    ]


def get_shortest_directed_path_between_subgraphs(graph: BELGraph, a: BELGraph, b: BELGraph) -> List[List[Any]]:
    """Calculate the shortest path(s) between disconnected sub-graphs ``a`` and ``b`` through ``graph``.

    :param graph: A BEL graph
    :param a: A sub-graph of :code:`graph`, disjoint from :code:`b`
    :param b: A sub-graph of :code:`graph`, disjoint from :code:`a`
    :return: A list of the shortest paths between the two sub-graphs
    """
    return _get_shortest_path_between_subgraphs_helper(graph, a, b)


def get_shortest_undirected_path_between_subgraphs(graph: BELGraph, a: BELGraph, b: BELGraph) -> List[List[Any]]:
    """Calculate the undirected shortest path(s) between disconnected sub-graphs ``a`` and ``b`` through ``graph``.

    :param graph: A BEL graph
    :param a: A sub-graph of :code:`graph`, disjoint from :code:`b`
    :param b: A sub-graph of :code:`graph`, disjoint from :code:`a`
    :return: A list of the shortest paths between the two sub-graphs
    """
    ug = graph.to_undirected()
    return _get_shortest_path_between_subgraphs_helper(ug, a, b)


def find_root_in_path(graph: BELGraph, path_nodes: List[BaseEntity]) -> Tuple[BELGraph, BaseEntity]:
    """Find the root of the path.

    This is defined as the node with the lowest out degree, if multiple:
    the root is the one with the highest out degree among those with lowest out degree

    :param graph: A BEL Graph
    :param path_nodes: A list of nodes in their order in a path
    :return: A pair of the graph: graph of the path and the root node
    """
    path_graph = graph.subgraph(path_nodes)

    # node_in_degree_tuple: list of tuples with (node,in_degree_of_node) in ascending order
    in_degrees = sorted(path_graph.in_degree().items(), key=itemgetter(1))

    # In case all have the same in degree it needs to be reference before
    tied_root_index = 0

    # Get index where the min in_degree stops (in case they are duplicates)
    for i in range(0, (len(in_degrees) - 1)):
        if in_degrees[i][1] < in_degrees[i + 1][1]:
            tied_root_index = i
            break

    # If there are multiple nodes with minimum in_degree take the one with max out degree
    # (in case multiple have the same out degree pick one random)
    if tied_root_index != 0:
        # node_out_degree_tuple: ordered list of tuples with (node,in_degree_of_node) in descending order
        out_degrees = sorted(path_graph.out_degree().items(), key=itemgetter(1), reverse=True)
        root_tuple = max(out_degrees[:tied_root_index], key=itemgetter(1))
    else:
        root_tuple = in_degrees[0]

    return path_graph, root_tuple[0]
