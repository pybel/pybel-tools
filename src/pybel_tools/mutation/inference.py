# -*- coding: utf-8 -*-

import logging

import itertools as itt

from pybel.constants import RELATION, TWO_WAY_RELATIONS
from pybel.struct import enrich_protein_and_rna_origins
from pybel.struct.pipeline import in_place_transformation, uni_in_place_transformation

__all__ = [
    'enrich_protein_and_rna_origins',
    'infer_missing_two_way_edges',
    'infer_missing_backwards_edge',
    'enrich_internal_unqualified_edges',
]

log = logging.getLogger(__name__)


@in_place_transformation
def infer_missing_two_way_edges(graph):
    """Add edges to the graph when a two way edge exists, and the opposite direction doesn't exist.

    Use: two way edges from BEL definition and/or axiomatic inverses of membership relations

    :param pybel.BELGraph graph: A BEL graph
    """
    for u, v, k, d in graph.edges(data=True, keys=True):
        if d[RELATION] in TWO_WAY_RELATIONS:
            infer_missing_backwards_edge(graph, u, v, k)


def infer_missing_backwards_edge(graph, u, v, k):
    """Add the same edge, but in the opposite direction if not already present.

    :type graph: pybel.BELGraph
    :type u: tuple
    :type v: tuple
    :type k: int
    """
    if u in graph[v]:
        for attr_dict in graph[v][u].values():
            if attr_dict == graph[u][v][k]:
                return

    graph.add_edge(v, u, key=k, **graph[u][v][k])


@uni_in_place_transformation
def enrich_internal_unqualified_edges(graph, subgraph):
    """Add the missing unqualified edges between entities in the subgraph that are contained within the full graph.

    :param pybel.BELGraph graph: The full BEL graph
    :param pybel.BELGraph subgraph: The query BEL subgraph
    """
    for u, v in itt.combinations(subgraph, 2):
        if not graph.has_edge(u, v):
            continue

        for k in graph[u][v]:
            if k < 0:
                subgraph.add_edge(u, v, key=k, **graph[u][v][k])
