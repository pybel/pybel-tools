# -*- coding: utf-8 -*-

import logging

import networkx as nx

from pybel import BELGraph
from pybel.struct.filters import filter_nodes, is_causal_relation
from pybel.struct.mutation import (
    get_multi_causal_downstream, get_multi_causal_upstream, get_random_subgraph, get_subgraph_by_all_shortest_paths,
    get_subgraph_by_annotation_value, get_subgraph_by_annotations, get_subgraph_by_authors, get_subgraph_by_edge_filter,
    get_subgraph_by_induction, get_subgraph_by_neighborhood, get_subgraph_by_pubmed, get_subgraph_by_second_neighbors,
)
from pybel.struct.operations import subgraph
from pybel.struct.pipeline import transformation
from pybel.struct.query import get_subgraph
from .search import search_node_names

log = logging.getLogger(__name__)

__all__ = [
    'get_subgraph_by_induction',
    'get_subgraph_by_node_filter',
    'get_subgraph_by_neighborhood',
    'get_subgraph_by_second_neighbors',
    'get_subgraph_by_all_shortest_paths',
    'get_subgraph_by_annotation_value',
    'get_subgraph_by_annotations',
    'get_subgraph_by_pubmed',
    'get_subgraph_by_authors',
    'get_subgraph_by_node_search',
    'get_causal_subgraph',
    'get_subgraph',
    'get_multi_causal_upstream',
    'get_multi_causal_downstream',
    'get_random_subgraph',
]


@transformation
def get_subgraph_by_node_filter(graph, node_filters):
    """Induces a graph on the nodes that pass all filters

    :param pybel.BELGraph graph: A BEL graph
    :param node_filters: A node filter or list/tuple of node filters
    :type node_filters: types.FunctionType or iter[types.FunctionType]
    :return: A subgraph induced over the nodes passing the given filters
    :rtype: pybel.BELGraph
    """
    return get_subgraph_by_induction(graph, filter_nodes(graph, node_filters))


@transformation
def get_causal_subgraph(graph: BELGraph) -> BELGraph:
    """Builds a new subgraph induced over all edges that are causal

    :param graph: A BEL graph
    :return: The causal sub-graph of the original BEL graph
    """
    return get_subgraph_by_edge_filter(graph, is_causal_relation)


@transformation
def get_subgraph_by_node_search(graph, query):
    """Gets a subgraph induced over all nodes matching the query string

    :param pybel.BELGraph graph: A BEL Graph
    :param str or iter[str] query: A query string or iterable of query strings for node names
    :return: A subgraph induced over the original BEL graph
    :rtype: pybel.BELGraph

    Thinly wraps :func:`search_node_names` and :func:`get_subgraph_by_induction`.
    """
    nodes = search_node_names(graph, query)
    return get_subgraph_by_induction(graph, nodes)


@transformation
def get_largest_component(graph):
    """Gets the giant component of a subgraph

    :param pybel.BELGraph graph: A BEL Graph
    :return: The giant component of the graph
    :rtype: pybel.BELGraph
    """
    biggest_component_nodes = max(nx.weakly_connected_components(graph), key=len)
    return subgraph(graph, biggest_component_nodes)
