# -*- coding: utf-8 -*-

import logging
import networkx as nx

from pybel.struct.filters import filter_nodes, is_causal_relation
from pybel.struct.mutation import (
    expand_nodes_neighborhoods, get_multi_causal_downstream, get_multi_causal_upstream,
    get_subgraph_by_all_shortest_paths, get_subgraph_by_annotation_value, get_subgraph_by_annotations,
    get_subgraph_by_authors, get_subgraph_by_edge_filter, get_subgraph_by_induction, get_subgraph_by_neighborhood,
    get_subgraph_by_pubmed, get_subgraph_by_second_neighbors,
)
from pybel.struct.pipeline import transformation
from .random_subgraph import get_random_subgraph
from .search import search_node_names
from ..filters.edge_filters import build_edge_data_filter

log = logging.getLogger(__name__)

__all__ = [
    'get_subgraph_by_induction',
    'get_subgraph_by_node_filter',
    'get_subgraph_by_neighborhood',
    'get_subgraph_by_second_neighbors',
    'get_subgraph_by_all_shortest_paths',
    'get_subgraph_by_annotation_value',
    'get_subgraph_by_annotations',
    'get_subgraph_by_data',
    'get_subgraph_by_pubmed',
    'get_subgraph_by_authors',
    'get_subgraph_by_node_search',
    'get_causal_subgraph',
    'get_subgraph',
    'get_multi_causal_upstream',
    'get_multi_causal_downstream',
    'get_random_subgraph',
]

#: Induce a subgraph over the given nodes
SEED_TYPE_INDUCTION = 'induction'
#: Induce a subgraph over the given nodes and expand to their first neighbors
SEED_TYPE_NEIGHBORS = 'neighbors'
#: Induce a subgraph over the given nodes and expand to their second neighbors
SEED_TYPE_DOUBLE_NEIGHBORS = 'dneighbors'
#: Induce a subgraph over the nodes in all shortest paths between the given nodes
SEED_TYPE_PATHS = 'shortest_paths'
#: Induce a subgraph over the edges provided by the given authors and their neighboring nodes
SEED_TYPE_AUTHOR = 'authors'
#: Induce a subgraph over the edges provided by the given citations and their neighboring nodes
SEED_TYPE_PUBMED = 'pubmed'
#: Generate an upstream candidate mechanism
SEED_TYPE_UPSTREAM = 'upstream'
#: Generate a downstream candidate mechanism
SEED_TYPE_DOWNSTREAM = 'downstream'
#: Induce a subgraph over the edges matching the given annotations
SEED_TYPE_ANNOTATION = 'annotation'
#: Induce a subgraph over a random set of (hopefully) connected edges
SEED_TYPE_SAMPLE = 'sample'

#: A set of the allowed seed type strings, as defined above
SEED_TYPES = {
    SEED_TYPE_INDUCTION,
    SEED_TYPE_NEIGHBORS,
    SEED_TYPE_DOUBLE_NEIGHBORS,
    SEED_TYPE_PATHS,
    SEED_TYPE_UPSTREAM,
    SEED_TYPE_DOWNSTREAM,
    SEED_TYPE_PUBMED,
    SEED_TYPE_AUTHOR,
    SEED_TYPE_ANNOTATION,
    SEED_TYPE_SAMPLE
}

#: Seed types that don't take node lists as their arguments
NONNODE_SEED_TYPES = {
    SEED_TYPE_ANNOTATION,
    SEED_TYPE_AUTHOR,
    SEED_TYPE_PUBMED,
    SEED_TYPE_SAMPLE,
}


class NodeDegreeIterError(ValueError):
    """Raised when failing to iterate over node degrees"""


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
def get_subgraph_by_data(graph, annotations):
    """Returns the subgraph filtering for Citation, Evidence or Annotation in the edges.
    
    :param pybel.BELGraph graph: A BEL graph
    :param dict annotations: Annotation filters (match all with :func:`pybel.utils.subdict_matches`)
    :return: A subgraph of the original BEL graph
    :rtype: pybel.BELGraph
    """
    return get_subgraph_by_edge_filter(graph, build_edge_data_filter(annotations))


@transformation
def get_causal_subgraph(graph):
    """Builds a new subgraph induced over all edges that are causal

    :param pybel.BELGraph graph: A BEL graph
    :return: A subgraph of the original BEL graph
    :rtype: pybel.BELGraph
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
def get_subgraph(graph, seed_method=None, seed_data=None, expand_nodes=None, remove_nodes=None):
    """Runs pipeline query on graph with multiple subgraph filters and expanders.

    Order of Operations:
    
    1. Seeding by given function name and data
    2. Add nodes
    3. Remove nodes

    :param pybel.BELGraph graph: A BEL graph
    :param str seed_method: The name of the get_subgraph_by_* function to use
    :param seed_data: The argument to pass to the get_subgraph function
    :param list[tuple] expand_nodes: Add the neighborhoods around all of these nodes
    :param list[tuple] remove_nodes: Remove these nodes and all of their in/out edges
    :rtype: Optional[pybel.BELGraph]
    """

    # Seed by the given function
    if seed_method == SEED_TYPE_INDUCTION:
        result = get_subgraph_by_induction(graph, seed_data)

    elif seed_method == SEED_TYPE_PATHS:
        result = get_subgraph_by_all_shortest_paths(graph, seed_data)

    elif seed_method == SEED_TYPE_NEIGHBORS:
        result = get_subgraph_by_neighborhood(graph, seed_data)

    elif seed_method == SEED_TYPE_DOUBLE_NEIGHBORS:
        result = get_subgraph_by_second_neighbors(graph, seed_data)

    elif seed_method == SEED_TYPE_UPSTREAM:
        result = get_multi_causal_upstream(graph, seed_data)

    elif seed_method == SEED_TYPE_DOWNSTREAM:
        result = get_multi_causal_downstream(graph, seed_data)

    elif seed_method == SEED_TYPE_PUBMED:
        result = get_subgraph_by_pubmed(graph, seed_data)

    elif seed_method == SEED_TYPE_AUTHOR:
        result = get_subgraph_by_authors(graph, seed_data)

    elif seed_method == SEED_TYPE_ANNOTATION:
        result = get_subgraph_by_annotations(graph, seed_data['annotations'], or_=seed_data.get('or'))

    elif seed_method == SEED_TYPE_SAMPLE:
        result = get_random_subgraph(
            graph,
            number_edges=seed_data.get('number_edges'),
            seed=seed_data.get('seed')
        )

    elif not seed_method:  # Otherwise, don't seed a subgraph
        result = graph.copy()
        log.debug('no seed function - using full network: %s', result.name)

    else:
        raise ValueError('Invalid seed method: {}'.format(seed_method))

    if result is None:
        log.debug('query returned no results')
        return

    log.debug('original graph has (%s nodes / %s edges)', result.number_of_nodes(), result.number_of_edges())

    # Expand around the given nodes
    if expand_nodes:
        expand_nodes_neighborhoods(graph, result, expand_nodes)
        log.debug('graph expanded to (%s nodes / %s edges)', result.number_of_nodes(), result.number_of_edges())

    # Delete the given nodes
    if remove_nodes:
        for node in remove_nodes:
            if node not in result:
                log.debug('%s is not in graph %s', node, graph.name)
                continue
            result.remove_node(node)
        log.debug('graph contracted to (%s nodes / %s edges)', result.number_of_nodes(), result.number_of_edges())

    log.debug(
        'Subgraph coming from %s (seed type) %s (data) contains %d nodes and %d edges',
        seed_method,
        seed_data,
        result.number_of_nodes(),
        result.number_of_edges()
    )

    return result


@transformation
def get_largest_component(graph):
    """Gets the giant component of a subgraph

    :param pybel.BELGraph graph: A BEL Graph
    :return: The giant component of the graph
    :rtype: pybel.BELGraph
    """
    return max(nx.weakly_connected_component_subgraphs(graph), key=len)
