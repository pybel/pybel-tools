# -*- coding: utf-8 -*-

"""This module contains functions that calculate properties of nodes"""

from collections import Counter

import networkx as nx

from pybel.constants import *
from pybel.struct.filters import get_nodes
from ..filters.node_filters import node_has_molecular_activity, node_is_degraded, node_is_translocated

__all__ = [
    'is_causal_relation',
    'get_causal_out_edges',
    'get_causal_in_edges',
    'has_causal_out_edges',
    'has_causal_in_edges',
    'is_causal_source',
    'is_causal_central',
    'is_causal_sink',
    'get_causal_source_nodes',
    'get_causal_central_nodes',
    'get_causal_sink_nodes',
    'get_degradations',
    'get_activities',
    'get_translocated',
    'count_variants',
    'count_top_centrality',
    'count_top_degrees',
]


def is_causal_relation(graph, u, v, k, d):
    """Is the given relation causal?"""
    return d[RELATION] in CAUSAL_RELATIONS


def get_causal_out_edges(graph, node):
    """Gets the out-edges to the given node that are causal

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :return: A set of (source, target) pairs where the source is the given node
    :rtype: set[tuple]
    """
    return {
        (u, v)
        for u, v, k, d in graph.out_edges_iter(node, keys=True, data=True)
        if is_causal_relation(graph, u, v, k, d)
    }


def get_causal_in_edges(graph, node):
    """Gets the in-edges to the given node that are causal

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :return: A set of (source, target) pairs where the target is the given node
    :rtype: set
    """
    return {
        (u, v)
        for u, v, k, d in graph.in_edges_iter(node, keys=True, data=True)
        if is_causal_relation(graph, u, v, k, d)
    }


def has_causal_out_edges(graph, node):
    """Does the node have causal out edges?

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :return: If the node has causal out edges
    :rtype: bool
    """
    return any(
        d[RELATION] in CAUSAL_RELATIONS
        for u, v, d in graph.out_edges_iter(node, data=True)
    )


def has_causal_in_edges(graph, node):
    """Does the node have causal in edges?

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :return: If the node has causal in edges
    :rtype: bool
    """
    return any(
        d[RELATION] in CAUSAL_RELATIONS
        for _, _, d in graph.in_edges_iter(node, data=True)
    )


def is_causal_source(graph, node):
    """Is the node is a causal source?

    - Doesn't have any causal in edge(s)
    - Does have causal out edge(s)

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :return: If the node is a causal source
    :rtype: bool
    """
    return not has_causal_in_edges(graph, node) and has_causal_out_edges(graph, node)


def is_causal_sink(graph, node):
    """Is the node is a causal sink?

    - Does have causal in edge(s)
    - Doesn't have any causal out edge(s)

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :return: If the node is a causal source
    :rtype: bool
    """
    return has_causal_in_edges(graph, node) and not has_causal_out_edges(graph, node)


def is_causal_central(graph, node):
    """Is the node neither a causal sink nor a causal source?

    - Does have causal in edges(s)
    - Does have causal out edge(s)

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :return: If the node is neither a causal sink nor a causal source
    :rtype: bool
    """
    return has_causal_in_edges(graph, node) and has_causal_out_edges(graph, node)


def get_causal_source_nodes(graph, function):
    """Returns a set of all nodes that have an in-degree of 0, which likely means that it is an external
    perturbagen and is not known to have any causal origin from within the biological system.

    These nodes are useful to identify because they generally don't provide any mechanistic insight.

    :param pybel.BELGraph graph: A BEL graph
    :param str function: The BEL function to filter by
    :return: A set of source nodes
    :rtype: set[tuple]
    """
    return {
        node
        for node, data in graph.nodes_iter(data=True)
        if data[FUNCTION] == function and is_causal_source(graph, node)
    }


def get_causal_central_nodes(graph, function):
    """Returns a set of all nodes that have both an in-degree > 0 and out-degree > 0. This means
    that they are an integral part of a pathway, since they are both produced and consumed.

    :param pybel.BELGraph graph: A BEL graph
    :param str function: The BEL function to filter by
    :return: A set of central ABUNDANCE nodes
    :rtype: set
    """
    return {
        node
        for node, data in graph.nodes_iter(data=True)
        if data[FUNCTION] == function and is_causal_central(graph, node)
    }


def get_causal_sink_nodes(graph, function):
    """Returns a set of all ABUNDANCE nodes that have an causal out-degree of 0, which likely means that the knowledge
    assembly is incomplete, or there is a curation error.

    :param pybel.BELGraph graph: A BEL graph
    :param str function: The BEL function to filter by
    :return: A set of sink ABUNDANCE nodes
    :rtype: set[tuple]
    """
    return {
        node
        for node, data in graph.nodes_iter(data=True)
        if data[FUNCTION] == function and is_causal_sink(graph, node)
    }


def get_degradations(graph):
    """Gets all nodes that are degraded

    :param pybel.BELGraph graph: A BEL graph
    :return: A set of nodes that are degraded
    :rtype: set[tuple]
    """
    return get_nodes(graph, node_is_degraded)


def get_activities(graph):
    """Gets all nodes that have molecular activities

    :param pybel.BELGraph graph: A BEL graph
    :return: A set of nodes that have molecular activities
    :rtype: set[tuple]
    """
    return get_nodes(graph, node_has_molecular_activity)


def get_translocated(graph):
    """Gets all nodes that are translocated

    :param pybel.BELGraph graph: A BEL graph
    :return: A set of nodes that are translocated
    :rtype: set[tuple]
    """
    return get_nodes(graph, node_is_translocated)


def node_has_variant(graph, node):
    """Does the node have variant information?

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :return: If the node contains variant information
    :rtype: bool
    """
    return VARIANTS in graph.node[node]


def count_variants(graph):
    """Counts how many of each type of variant a graph has"""
    return Counter(
        variant_data[KIND]
        for node, data in graph.nodes_iter(data=True)
        if node_has_variant(graph, node)
        for variant_data in data[VARIANTS]
    )


def count_top_degrees(graph, number=30):
    dd = graph.degree()
    dc = Counter(dd)
    return dict(dc.most_common(number))


def count_top_centrality(graph, number=30):
    dd = nx.betweenness_centrality(graph)
    dc = Counter(dd)
    return dict(dc.most_common(number))
