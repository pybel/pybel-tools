# -*- coding: utf-8 -*-

import logging

from pybel.constants import (
    EQUIVALENT_TO, FUNCTION, GENE, HAS_VARIANT, ORTHOLOGOUS, PROTEIN, RELATION, TRANSCRIBED_TO,
    VARIANTS,
)
from pybel.struct.filters import (
    filter_edges, has_polarity,
)
from pybel.struct.mutation import (
    collapse_nodes, collapse_pair, collapse_to_genes, get_subgraph_by_edge_filter,
)
from pybel.struct.pipeline import in_place_transformation, transformation
from ..filters.edge_filters import build_relation_filter, build_source_namespace_filter, build_target_namespace_filter
from ..summary.edge_summary import pair_is_consistent

__all__ = [
    'collapse_nodes',
    'rewire_variants_to_genes',
    'collapse_all_variants',
    'collapse_gene_variants',
    'collapse_protein_variants',
    'collapse_consistent_edges',
    'collapse_equivalencies_by_namespace',
    'collapse_orthologies_by_namespace',
    'collapse_to_protein_interactions',
]

log = logging.getLogger(__name__)


def _collapse_variants_by_function(graph, func):
    """Collapses all of the given functions' variants' edges to their parents, in-place

    :param pybel.BELGraph graph: A BEL graph
    :param str func: A BEL function
    """
    for parent_node, variant_node, data in graph.edges(data=True):
        if data[RELATION] == HAS_VARIANT and graph.node[parent_node][FUNCTION] == func:
            collapse_pair(graph, from_node=variant_node, to_node=parent_node)


@in_place_transformation
def collapse_protein_variants(graph):
    """Collapses all protein's variants' edges to their parents, in-place

    :param pybel.BELGraph graph: A BEL graph
    """
    _collapse_variants_by_function(graph, PROTEIN)


@in_place_transformation
def collapse_gene_variants(graph):
    """Collapses all gene's variants' edges to their parents, in-place

    :param pybel.BELGraph graph: A BEL graph
    """
    _collapse_variants_by_function(graph, GENE)


@in_place_transformation
def rewire_variants_to_genes(graph):
    """Finds all protein variants that are pointing to a gene and not a protein and fixes them by changing their
    function to be :data:`pybel.constants.GENE`, in place

    :param pybel.BELGraph graph: A BEL graph
    
    A use case is after running :func:`collapse_to_genes`.
    """
    for node, data in graph.nodes(data=True):
        if data[FUNCTION] != PROTEIN:
            continue
        if VARIANTS not in data:
            continue
        if any(d[RELATION] == TRANSCRIBED_TO for u, v, d in graph.in_edges_iter(data=True)):
            graph.node[node][FUNCTION] = GENE


@in_place_transformation
def collapse_all_variants(graph):
    """Collapse all genes', RNAs', miRNAs', and proteins' variants to their parents.
    
    :param pybel.BELGraph graph: A BEL Graph
    """
    for parent_node, variant_node, d in graph.edges(data=True):
        if d[RELATION] == HAS_VARIANT:
            collapse_pair(graph, survivor=parent_node, victim=variant_node)


def _collapse_edge_passing_predicates(graph, edge_predicates=None):
    """Collapses all edges passing the given edge predicates

    :param pybel.BELGraph graph: A BEL Graph
    :param edge_predicates: A predicate or list of predicates
    :type edge_predicates: Optional[(pybel.BELGraph, tuple, tuple, int) -> bool or iter[(pybel.BELGraph, tuple, tuple, int) -> bool]]
    """
    for u, v, _ in filter_edges(graph, edge_predicates=edge_predicates):
        collapse_pair(graph, survivor=u, victim=v)


def _collapse_edge_by_namespace(graph, victim_namespace, survivor_namespace, relation):
    """Collapses pairs of nodes with the given namespaces that have the given relationship

    :param pybel.BELGraph graph: A BEL Graph
    :param str or iter[str] victim_namespace: The namespace(s) of the node to collapse
    :param str or survivor_namespace: The namespace of the node to keep
    :param relation: The relation to search
    :type relation: str or iter[str]
    """
    relation_filter = build_relation_filter(relation)
    source_namespace_filter = build_source_namespace_filter(victim_namespace)
    target_namespace_filter = build_target_namespace_filter(survivor_namespace)

    edge_predicates = [
        relation_filter,
        source_namespace_filter,
        target_namespace_filter
    ]

    _collapse_edge_passing_predicates(graph, edge_predicates=edge_predicates)


@in_place_transformation
def collapse_equivalencies_by_namespace(graph, victim_namespace, survivor_namespace):
    """Collapse pairs of nodes with the given namespaces that have equivalence relationships.
    
    :param pybel.BELGraph graph: A BEL graph
    :param str or iter[str] victim_namespace: The namespace(s) of the node to collapse
    :param str survivor_namespace: The namespace of the node to keep

    To convert all ChEBI names to InChI keys, assuming there are appropriate equivalence relations between nodes with
    those namespaces:
    
    >>> collapse_equivalencies_by_namespace(graph, 'CHEBI', 'CHEBIID')
    >>> collapse_equivalencies_by_namespace(graph, 'CHEBIID', 'INCHI')
    """
    _collapse_edge_by_namespace(graph, victim_namespace, survivor_namespace, EQUIVALENT_TO)


@in_place_transformation
def collapse_orthologies_by_namespace(graph, victim_namespace, survivor_namespace):
    """Collapse pairs of nodes with the given namespaces that have orthology relationships.

    :param pybel.BELGraph graph: A BEL Graph
    :param str or iter[str] victim_namespace: The namespace(s) of the node to collapse
    :param str survivor_namespace: The namespace of the node to keep

    To collapse all MGI nodes to their HGNC orthologs, use:
    >>> collapse_orthologies_by_namespace('MGI', 'HGNC')


    To collapse collapse both MGI and RGD nodes to their HGNC orthologs, use:
    >>> collapse_orthologies_by_namespace(['MGI', 'RGD'], 'HGNC')
    """
    _collapse_edge_by_namespace(graph, victim_namespace, survivor_namespace, ORTHOLOGOUS)


@in_place_transformation
def collapse_entrez_to_hgnc(graph):
    """Collapse Entrez equivalences to HGNC.

    :param pybel.BELGraph graph: A BEL graph
    """
    collapse_equivalencies_by_namespace(graph, ['EGID', 'EG', 'ENTREZ'], 'HGNC')


@in_place_transformation
def collapse_mgi_to_hgnc(graph):
    """Collapse MGI orthologies to HGNC.

    :param pybel.BELGraph graph: A BEL graph
    """
    collapse_orthologies_by_namespace(graph, ['MGI', 'MGIID'], 'HGNC')


@in_place_transformation
def collapse_rgd_to_hgnc(graph):
    """Collapse RGD orthologies to HGNC.

    :param pybel.BELGraph graph: A BEL graph
    """
    collapse_orthologies_by_namespace(graph, ['RGD', 'RGDID'], 'HGNC')


@in_place_transformation
def collapse_flybase_to_hgnc(graph):
    """Collapse FlyBase orthologies to HGNC.

    :param pybel.BELGraph graph: A BEL graph
    """
    collapse_orthologies_by_namespace(graph, 'FLYBASE', 'HGNC')


@in_place_transformation
def collapse_entrez_equivalencies(graph):
    """Collapses all equivalence edges away from Entrez. Assumes well formed, 2-way equivalencies

    :param pybel.BELGraph graph: A BEL graph
    """
    relation_filter = build_relation_filter(EQUIVALENT_TO)
    source_namespace_filter = build_source_namespace_filter(['EGID', 'EG', 'ENTREZ'])

    edge_predicates = [
        relation_filter,
        source_namespace_filter,
    ]

    _collapse_edge_passing_predicates(graph, edge_predicates=edge_predicates)


@in_place_transformation
def collapse_consistent_edges(graph):
    """Collapse consistent edges together.

    .. warning:: This operation doesn't preserve evidences or other annotations

    :param pybel.BELGraph graph: A BEL Graph
    """
    for u, v in graph.edges():
        relation = pair_is_consistent(graph, u, v)

        if not relation:
            continue

        edges = [(u, v, k) for k in graph[u][v]]
        graph.remove_edges_from(edges)
        graph.add_edge(u, v, attr_dict={RELATION: relation})


@transformation
def collapse_to_protein_interactions(graph):
    """Collapse to a graph made of only causal gene/protein edges.

    :param pybel.BELGraph graph: A BEL Graph
    """
    rv = graph.copy()

    collapse_to_genes(rv)

    def is_edge_ppi(g, u, v, k):
        return g.node[u][FUNCTION] == GENE and g.node[v][FUNCTION] == GENE

    return get_subgraph_by_edge_filter(rv, edge_predicates=[has_polarity, is_edge_ppi])
