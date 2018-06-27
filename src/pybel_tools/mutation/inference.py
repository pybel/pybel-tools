# -*- coding: utf-8 -*-

import itertools as itt

import logging

from pybel.constants import RELATION, TWO_WAY_RELATIONS, unqualified_edge_code
from pybel.struct import enrich_protein_and_rna_origins
from pybel.struct.filters import filter_edges
from pybel.struct.pipeline import in_place_transformation, uni_in_place_transformation
from ..constants import INFERRED_INVERSE
from ..filters.edge_filters import build_relation_filter
from ..utils import safe_add_edge

__all__ = [
    'enrich_protein_and_rna_origins',
    'infer_missing_two_way_edges',
    'infer_missing_backwards_edge',
    'infer_missing_inverse_edge',
    'enrich_internal_unqualified_edges',
]

log = logging.getLogger(__name__)


@in_place_transformation
def infer_missing_two_way_edges(graph):
    """If a two way edge exists, and the opposite direction doesn't exist, add it to the graph

    Use: two way edges from BEL definition and/or axiomatic inverses of membership relations

    :param pybel.BELGraph graph: A BEL graph
    """
    for u, v, k, d in graph.edges_iter(data=True, keys=True):
        if d[RELATION] in TWO_WAY_RELATIONS:
            infer_missing_backwards_edge(graph, u, v, k)


@in_place_transformation
def infer_missing_inverse_edge(graph, relations):
    """Adds inferred edges based on pre-defined axioms

    :param pybel.BELGraph graph: A BEL network
    :param relations: single or iterable of relation names to add their inverse inferred edges
    :type relations: str or iter[str]
    """

    if isinstance(relations, str):
        return infer_missing_inverse_edge(graph, [relations])

    for u, v, _, d in filter_edges(graph, build_relation_filter(relations)):
        relation = d[RELATION]
        graph.add_edge(v, u, key=unqualified_edge_code[relation], **{RELATION: INFERRED_INVERSE[relation]})


@in_place_transformation
def infer_missing_backwards_edge(graph, u, v, k):
    """Adds the same edge, but in the opposite direction if not already present

    :param pybel.BELGraph graph: A BEL graph
    :param u: A BEL node
    :type u: tuple
    :param v: A BEL node
    :type v: tuple
    :param k: The edge key
    :type k: int
    """
    if u in graph.edge[v]:
        for attr_dict in graph.edge[v][u].values():
            if attr_dict == graph.edge[u][v][k]:
                return

    safe_add_edge(graph, v, u, key=k, attr_dict=graph.edge[u][v][k])


@uni_in_place_transformation
def enrich_internal_unqualified_edges(graph, subgraph):
    """Adds the missing unqualified edges between entities in the subgraph that are contained within the full graph

    :param pybel.BELGraph graph: The full BEL graph
    :param pybel.BELGraph subgraph: The query BEL subgraph
    """
    for u, v in itt.combinations(subgraph, 2):
        if not graph.has_edge(u, v):
            continue

        for k in graph.edge[u][v]:
            if k < 0:
                subgraph.add_edge(u, v, key=k, attr_dict=graph.edge[u][v][k])
