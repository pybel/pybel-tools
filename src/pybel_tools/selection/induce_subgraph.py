# -*- coding: utf-8 -*-

import networkx as nx

from pybel import BELGraph
from pybel.struct.filters import filter_nodes, is_causal_relation
from pybel.struct.filters.typing import NodePredicates
from pybel.struct.mutation import get_subgraph_by_edge_filter, get_subgraph_by_induction
from pybel.struct.operations import subgraph
from pybel.struct.pipeline import transformation
from pybel.typing import Strings
from .search import search_node_names

__all__ = [
    'get_subgraph_by_node_filter',
    'get_causal_subgraph',
    'get_subgraph_by_node_search',
    'get_largest_component',
]


@transformation
def get_subgraph_by_node_filter(graph: BELGraph, node_predicates: NodePredicates) -> BELGraph:
    """Induce a sub-graph on the nodes that pass the given predicate(s)."""
    return get_subgraph_by_induction(graph, filter_nodes(graph, node_predicates))


@transformation
def get_causal_subgraph(graph: BELGraph) -> BELGraph:
    """Build a new sub-graph induced over the causal edges."""
    return get_subgraph_by_edge_filter(graph, is_causal_relation)


@transformation
def get_subgraph_by_node_search(graph: BELGraph, query: Strings) -> BELGraph:
    """Get a sub-graph induced over all nodes matching the query string.

    :param graph: A BEL Graph
    :param query: A query string or iterable of query strings for node names

    Thinly wraps :func:`search_node_names` and :func:`get_subgraph_by_induction`.
    """
    nodes = search_node_names(graph, query)
    return get_subgraph_by_induction(graph, nodes)


@transformation
def get_largest_component(graph: BELGraph) -> BELGraph:
    """Get the giant component of a graph."""
    biggest_component_nodes = max(nx.weakly_connected_components(graph), key=len)
    return subgraph(graph, biggest_component_nodes)
