# -*- coding: utf-8 -*-

from .. import pipeline

__all__ = [
    'NODE_HIGHLIGHT',
    'EDGE_HIGHLIGHT',
    'is_node_highlighted',
    'highlight_nodes',
    'remove_highlight_nodes',
    'is_edge_highlighted',
    'highlight_edges',
    'remove_highlight_edges',
    'highlight_subgraph',
    'remove_highlight_subgraph',
]

NODE_HIGHLIGHT = 'pybel_highlight'
NODE_HIGHLIGHT_DEFAULT_COLOR = 'orange'
EDGE_HIGHLIGHT = 'pybel_highlight'
EDGE_HIGHLIGHT_DEFAULT_COLOR = 'orange'


@pipeline.in_place_mutator
def highlight_nodes(graph, nodes=None, color=None):
    """Adds a highlight tag to the given nodes.

    :param pybel.BELGraph graph: A BEL graph
    :param nodes: The nodes to add a highlight tag on
    :type nodes: iter[tuple]
    :param color: The color to highlight (use something that works with CSS)
    :type color: str
    """
    for node in nodes if nodes is not None else graph.nodes_iter():
        graph.node[node][NODE_HIGHLIGHT] = NODE_HIGHLIGHT_DEFAULT_COLOR if color is None else color


def is_node_highlighted(graph, node):
    """Returns if the given node is highlighted.

    :param pybel.BELGraph graph: A BEL graph
    :param node: A BEL node
    :type node: tuple
    :return: Does the node contain highlight information?
    :rtype: bool
    """
    return NODE_HIGHLIGHT in graph.node[node]


@pipeline.in_place_mutator
def remove_highlight_nodes(graph, nodes=None):
    """Removes the highlight from the given nodes, or all nodes if none given.

    :param pybel.BELGraph graph: A BEL graph
    :param nodes: The list of nodes to un-highlight
    :type nodes: list
    """
    for node in graph.nodes_iter() if nodes is None else nodes:
        if is_node_highlighted(graph, node):
            del graph.node[node][NODE_HIGHLIGHT]


@pipeline.in_place_mutator
def highlight_edges(graph, edges=None, color=None):
    """Adds a highlight tag to the given edges.

    :param pybel.BELGraph graph: A BEL graph
    :param edges: The edges (4-tuples of u, v, k, d) to add a highlight tag on
    :type edges: iter[tuple]
    :param str color: The color to highlight (use something that works with CSS)
    """
    for u, v, k, d in edges if edges is not None else graph.edges_iter(keys=True, data=True):
        graph.edge[u][v][k][EDGE_HIGHLIGHT] = EDGE_HIGHLIGHT_DEFAULT_COLOR if color is None else color


def is_edge_highlighted(graph, u, v, k, d):
    """Returns if the given edge is highlighted.
    
    :param pybel.BELGraph graph: A BEL graph
    :return: Does the edge contain highlight information?
    :rtype: bool
    """
    return EDGE_HIGHLIGHT in graph.edge[u][v][k]


@pipeline.in_place_mutator
def remove_highlight_edges(graph, edges=None):
    """Removes the highlight from the given edges, or all edges if none given.

    :param pybel.BELGraph graph: A BEL graph
    :param edges: The edges (4-tuple of u,v,k,d) to remove the highlight from)
    :type edges: iter[tuple]
    """
    for u, v, k, d in graph.edges_iter(keys=True, data=True) if edges is None else edges:
        if is_edge_highlighted(graph, u, v, k, d):
            del graph.edge[u][v][k][EDGE_HIGHLIGHT]


@pipeline.uni_in_place_mutator
def highlight_subgraph(universe, graph):
    """Highlights all nodes/edges in the universe that in the given graph.

    :param pybel.BELGraph universe: The universe of knowledge
    :type pybel.BELGraph graph: The BEL graph to mutate
    """
    highlight_nodes(universe, graph.nodes_iter())
    highlight_edges(universe, graph.edges_iter())


@pipeline.in_place_mutator
def remove_highlight_subgraph(graph, subgraph):
    """Removes the highlight from all nodes/edges in the graph that are in the subgraph.

    :param pybel.BELGraph graph: The BEL graph to mutate
    :param pybel.BELGraph subgraph: The subgraph from which to remove the highlighting
    """
    remove_highlight_nodes(graph, subgraph.nodes_iter())
    remove_highlight_edges(graph, subgraph.edges_iter())
