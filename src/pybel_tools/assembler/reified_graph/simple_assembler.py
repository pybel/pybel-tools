# -*- coding: utf-8 -*-

"""A simple assembler that doesn't do any interpreting in the middle."""

import networkx as nx

from pybel import BELGraph

__all__ = [
    'reify_bel_graph_simple',
]


def reify_bel_graph_simple(graph: BELGraph) -> nx.DiGraph:
    """Convert a BEL graph to a simple reified graph.

    In this graph, there are two types of nodes:

    1. BEL Terms
    2. BEL Edges

    There are two types of edges:

    1. The subject is a BEL term and the object is a BEL edge. This represents
       the subject of the BEL statement.
    2. The subject is a BEL edge and the object is a BEL term. This represents
       the object of the BEL statement.
    """
    rv = nx.DiGraph()
    for u, v, edge_data in graph.edges(data=True):
        bel = graph.edge_to_bel(u, v, edge_data)
        rv.add_edge(u, bel)
        rv.add_edge(bel, v)
    return rv
