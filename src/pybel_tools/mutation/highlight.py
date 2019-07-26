# -*- coding: utf-8 -*-

"""Functions for adding highlighting tags to nodes and edges in BEL graphs."""

from typing import Iterable, Optional

from pybel import BELGraph
from pybel.dsl import BaseEntity
from pybel.struct.filters.typing import EdgeIterator
from pybel.struct.pipeline import in_place_transformation, uni_in_place_transformation

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


@in_place_transformation
def highlight_nodes(
        graph: BELGraph,
        nodes: Optional[Iterable[BaseEntity]] = None,
        color: Optional[str] = None,
) -> None:
    """Add a highlight tag to the given nodes.

    :param graph: A BEL graph
    :param nodes: The nodes to add a highlight tag on
    :param color: The color to highlight (use something that works with CSS)
    """
    color = color or NODE_HIGHLIGHT_DEFAULT_COLOR
    for node in nodes if nodes is not None else graph:
        graph.nodes[node][NODE_HIGHLIGHT] = color


def is_node_highlighted(graph: BELGraph, node: BaseEntity) -> bool:
    """Return if the given node is highlighted."""
    return NODE_HIGHLIGHT in graph.nodes[node]


@in_place_transformation
def remove_highlight_nodes(
        graph: BELGraph,
        nodes: Optional[Iterable[BaseEntity]] = None,
) -> None:
    """Remove the highlight from the given nodes, or all nodes if none given.

    :param graph: A BEL graph
    :param nodes: The list of nodes to un-highlight
    """
    for node in graph if nodes is None else nodes:
        if is_node_highlighted(graph, node):
            del graph.nodes[node][NODE_HIGHLIGHT]


@in_place_transformation
def highlight_edges(
        graph: BELGraph,
        edges: Optional[EdgeIterator] = None,
        color: Optional[str] = None,
) -> None:
    """Add a highlight tag to the given edges.

    :param graph: A BEL graph
    :param edges: The edges (3-tuples of u, v, k) to add a highlight tag on
    :param color: The color to highlight (use something that works with CSS)
    """
    if color is None:
        color = EDGE_HIGHLIGHT_DEFAULT_COLOR

    for u, v, k, d in edges if edges is not None else graph.edges(keys=True, data=True):
        graph[u][v][k][EDGE_HIGHLIGHT] = color


def is_edge_highlighted(graph: BELGraph, u: BaseEntity, v: BaseEntity, k: str) -> bool:
    """Return if the given edge is highlighted."""
    return EDGE_HIGHLIGHT in graph[u][v][k]


@in_place_transformation
def remove_highlight_edges(
        graph: BELGraph,
        edges: Optional[EdgeIterator] = None,
) -> None:
    """Remove the highlight from the given edges, or all edges if none given.

    :param graph: A BEL graph
    :param edges: The edges (3-tuple of u, v, k) to remove the highlight from)
    """
    for u, v, k in graph.edges(keys=True) if edges is None else edges:
        if is_edge_highlighted(graph, u, v, k):
            del graph[u][v][k][EDGE_HIGHLIGHT]


@uni_in_place_transformation
def highlight_subgraph(universe: BELGraph, graph: BELGraph) -> None:
    """Highlight all nodes/edges in the universe that in the given graph.

    :param universe: The universe of knowledge
    :param graph: The BEL graph to mutate
    """
    highlight_nodes(universe, graph.nodes())
    highlight_edges(universe, graph.edges())


@in_place_transformation
def remove_highlight_subgraph(graph: BELGraph, subgraph: BELGraph) -> None:
    """Remove the highlight from all nodes/edges in the graph that are in the subgraph.

    :param graph: The BEL graph to mutate
    :param subgraph: The subgraph from which to remove the highlighting
    """
    remove_highlight_nodes(graph, subgraph.nodes())
    remove_highlight_edges(graph, subgraph.edges())
