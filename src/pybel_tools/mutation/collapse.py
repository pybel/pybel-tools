# -*- coding: utf-8 -*-

import logging

from pybel.constants import *
from pybel.struct.filters import get_nodes
from pybel.struct.filters.edge_filters import filter_edges, invert_edge_predicate
from pybel.struct.filters.edge_predicates import has_polarity
from pybel.struct.filters.node_predicate_builders import function_inclusion_filter_builder
from pybel.struct.mutation.collapse import (
    build_central_dogma_collapse_dict, build_central_dogma_collapse_gene_dict,
    collapse_by_central_dogma, collapse_by_central_dogma_to_genes, collapse_nodes, collapse_pair,
)
from pybel.struct.mutation.inference import infer_central_dogma
from pybel.struct.pipeline import in_place_transformation, transformation
from pybel.struct.utils import update_node_helper
from ..filters.edge_filters import build_relation_filter, build_source_namespace_filter, build_target_namespace_filter
from ..mutation.deletion import remove_filtered_edges
from ..summary.edge_summary import pair_is_consistent
from ..utils import all_edges_iter

__all__ = [
    'collapse_nodes',
    'build_central_dogma_collapse_dict',
    'build_central_dogma_collapse_gene_dict',
    'collapse_by_central_dogma_to_genes',
    'collapse_by_central_dogma_to_genes_out_place',
    'rewire_variants_to_genes',
    'collapse_all_variants',
    'collapse_all_variants_out_place',
    'collapse_gene_variants',
    'collapse_protein_variants',
    'collapse_consistent_edges',
    'collapse_equivalencies_by_namespace',
    'collapse_orthologies_by_namespace',
    'collapse_to_protein_interactions',
]

log = logging.getLogger(__name__)


@transformation
def collapse_by_central_dogma_to_genes_out_place(graph):
    """Collapses all nodes from the central dogma (:data:`pybel.constants.GENE`, :data:`pybel.constants.RNA`, 
    :data:`pybel.constants.MIRNA`, and :data:`pybel.constants.PROTEIN`) to :data:`pybel.constants.GENE`, in-place. This 
    function wraps :func:`collapse_nodes` and :func:`build_central_dogma_collapse_gene_dict`.

    :param pybel.BELGraph graph: A BEL graph
    :rtype: pybel.BELGraph
    
    Equivalent to:
    
    >>> infer_central_dogma(graph)
    >>> collapse_nodes(graph, build_central_dogma_collapse_gene_dict(graph))
    """
    result = graph.copy()
    collapse_by_central_dogma_to_genes(result)
    return result


@in_place_transformation
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
    
    A use case is after running :func:`collapse_by_central_dogma_to_genes`
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
    """Collapses all ``hasVariant`` edges to the parent node, in place
    
    :param pybel.BELGraph graph: A BEL Graph
    """
    for parent_node, variant_node, d in graph.edges(data=True):
        if d[RELATION] == HAS_VARIANT:
            collapse_pair(graph, from_node=variant_node, to_node=parent_node)


@transformation
def collapse_all_variants_out_place(graph):
    """Collapses all ``hasVariant`` edges to the parent node, not in place

    :param pybel.BELGraph graph: A BEL Graph
    :rtype: pybel.BELGraph
    """
    result = graph.copy()
    _collapse_variants_by_function(result)
    return result


@in_place_transformation
def collapse_by_opening_on_central_dogma(graph):
    """Infers the matching RNA for each protein and the gene for each RNA and miRNA, then collapses the corresponding
    gene node to its RNA/miRNA node, then possibly from RNA to protein if available. Wraps :func:`infer_central_dogma`
    and :func:`collapse_by_central_dogma`.

    :param pybel.BELGraph graph: A BEL graph

    Equivalent to:

    >>> infer_central_dogma(graph)
    >>> collapse_by_central_dogma(graph)
    """
    infer_central_dogma(graph)
    collapse_by_central_dogma(graph)


@in_place_transformation
def collapse_by_opening_by_central_dogma_to_genes(graph):
    """Infers the matching RNA for each protein and the gene for each RNA and miRNA, then collapses the corresponding
    protein and RNA/miRNA nodes to the gene node.

    This method is useful to help overcome issues with BEL curation, when curators sometimes haphazardly annotate
    entities as either a gene, RNA, or protein. There is possibly significant biological subtlty that can be lost
    during this process, but sometimes this must be done to overcome the noise introduced by these kinds of mistakes.
    
    Wraps :func:`infer_central_dogma` and :func:`collapse_by_central_dogma_to_genes`.
    
    :param pybel.BELGraph graph: A BEL graph

    Equivalent to:

    >>> infer_central_dogma(graph)
    >>> collapse_by_central_dogma_to_genes(graph)
    """
    infer_central_dogma(graph)
    collapse_by_central_dogma_to_genes(graph)


def _collapse_edge_passing_predicates(graph, edge_predicates=None):
    """Collapses all edges passing the given edge predicates

    :param pybel.BELGraph graph: A BEL Graph
    :param edge_predicates: A predicate or list of predicates
    :type edge_predicates: Optional[(pybel.BELGraph, tuple, tuple, int) -> bool or iter[(pybel.BELGraph, tuple, tuple, int) -> bool]]
    """
    for u, v, _ in filter_edges(graph, edge_predicates=edge_predicates):
        collapse_pair(graph, from_node=u, to_node=v)


def _collapse_edge_by_namespace(graph, from_namespace, to_namespace, relation):
    """Collapses pairs of nodes with the given namespaces that have the given relationship

    :param pybel.BELGraph graph: A BEL Graph
    :param str or iter[str] from_namespace: The namespace of the node to collapse
    :param str or iter[str] to_namespace: The namespace of the node to keep
    :param relation: The relation to search
    :type relation: str or iter[str]
    """
    relation_filter = build_relation_filter(relation)
    source_namespace_filter = build_source_namespace_filter(from_namespace)
    target_namespace_filter = build_target_namespace_filter(to_namespace)

    edge_predicates = [
        relation_filter,
        source_namespace_filter,
        target_namespace_filter
    ]

    _collapse_edge_passing_predicates(graph, edge_predicates=edge_predicates)


@in_place_transformation
def collapse_equivalencies_by_namespace(graph, from_namespace, to_namespace):
    """Collapses pairs of nodes with the given namespaces that have equivalence relationships
    
    :param pybel.BELGraph graph: A BEL graph
    :param str or iter[str] from_namespace: The namespace of the node to collapse
    :param str or iter[str] to_namespace: The namespace of the node to keep

    To convert all ChEBI names to InChI keys, assuming there are appropriate equivalence relations between nodes with
    those namespaces:
    
    >>> collapse_equivalencies_by_namespace(graph, 'CHEBI', 'CHEBIID')
    >>> collapse_equivalencies_by_namespace(graph, 'CHEBIID', 'INCHI')
    """
    _collapse_edge_by_namespace(graph, from_namespace, to_namespace, EQUIVALENT_TO)


@in_place_transformation
def collapse_orthologies_by_namespace(graph, from_namespace, to_namespace):
    """Collapses pairs of nodes with the given namespaces that have orthology relationships

    :param pybel.BELGraph graph: A BEL Graph
    :param str or iter[str] from_namespace: The namespace of the node to collapse
    :param str or iter[str] to_namespace: The namespace of the node to keep

    To collapse all MGI nodes to their HGNC orthologs, use:

    >>> collapse_orthologies_by_namespace('MGI', 'HGNC')
    """
    _collapse_edge_by_namespace(graph, from_namespace, to_namespace, ORTHOLOGOUS)


@in_place_transformation
def collapse_entrez_to_hgnc(graph):
    """Collapses Entrez orthologies to HGNC.

    Implemented with :func:`collapse_orthologies_by_namespace`.

    :param pybel.BELGraph graph: A BEL graph
    """
    collapse_orthologies_by_namespace(graph, ['EGID', 'EG', 'ENTREZ'], 'HGNC')


@in_place_transformation
def collapse_mgi_to_hgnc(graph):
    """Collapses MGI orthologies to HGNC

    Implemented with :func:`collapse_orthologies_by_namespace`.

    :param pybel.BELGraph graph: A BEL graph
    """
    collapse_orthologies_by_namespace(graph, ['MGI', 'MGIID'], 'HGNC')


@in_place_transformation
def collapse_rgd_to_hgnc(graph):
    """Collapses RGD orthologies to HGNC

    Implemented with :func:`collapse_orthologies_by_namespace`.

    :param pybel.BELGraph graph: A BEL graph
    """
    collapse_orthologies_by_namespace(graph, ['RDG', 'RGDID'], 'HGNC')


@in_place_transformation
def collapse_flybase_to_hgnc(graph):
    """Collapses FlyBase orthologies to HGNC

    Implemented with :func:`collapse_orthologies_by_namespace`.

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
    """Collapses consistent edges together

    .. warning:: This operation doesn't preserve evidences or other annotations

    :param pybel.BELGraph graph: A BEL Graph
    """
    for u, v in graph.edges():
        relation = pair_is_consistent(graph, u, v)

        if not relation:
            continue

        edges = list(all_edges_iter(graph, u, v))
        graph.remove_edges_from(edges)
        graph.add_edge(u, v, attr_dict={RELATION: relation})


@transformation
def collapse_to_protein_interactions(graph):
    """Collapse to a graph made of only causal gene/protein edges.

    :param pybel.BELGraph graph: A BEL Graph
    """
    collapse_by_central_dogma_to_genes(graph)

    remove_filtered_edges(graph, invert_edge_predicate(has_polarity))

    filtered_graph = graph.subgraph(get_nodes(graph, node_predicates=function_inclusion_filter_builder(GENE)))
    update_node_helper(graph, filtered_graph)

    return filtered_graph
