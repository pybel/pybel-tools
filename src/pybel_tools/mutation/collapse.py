# -*- coding: utf-8 -*-

import logging
from collections import defaultdict

from pybel.constants import *
from pybel.struct.filters import get_nodes
from pybel.struct.filters.edge_filters import filter_edges, invert_edge_filter
from pybel.struct.filters.edge_predicates import has_polarity
from .inference import infer_central_dogma
from .. import pipeline
from ..filters.edge_filters import build_relation_filter, build_source_namespace_filter, build_target_namespace_filter
from ..filters.node_filters import function_inclusion_filter_builder
from ..mutation.deletion import remove_filtered_edges
from ..summary.edge_summary import pair_is_consistent
from ..utils import all_edges_iter

__all__ = [
    'collapse_nodes',
    'build_central_dogma_collapse_dict',
    'build_central_dogma_collapse_gene_dict',
    'collapse_by_central_dogma',
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


@pipeline.in_place_mutator
def collapse_pair(graph, *, to_node, from_node):
    """Rewires all edges from the synonymous node to the survivor node, then deletes the synonymous node.
    
    Does not keep edges between the two nodes.
    
    :param pybel.BELGraph graph: A BEL graph
    :param tuple to_node: The BEL node to collapse all edges on the synonym to
    :param tuple from_node: The BEL node to collapse into the surviving node
    """
    for successor in graph.successors_iter(from_node):
        if successor == to_node:
            continue

        for key, data in graph.edge[from_node][successor].items():
            if key >= 0:  # FIXME switch to safe add edge?
                graph.add_edge(to_node, successor, attr_dict=data)
            elif successor not in graph.edge[to_node] or key not in graph.edge[to_node][successor]:
                graph.add_edge(to_node, successor, key=key, **{RELATION: unqualified_edges[-1 - key]})

    for predecessor in graph.predecessors_iter(from_node):
        if predecessor == to_node:
            continue

        for key, data in graph.edge[predecessor][from_node].items():
            if key >= 0:  # FIXME switch to safe add edge?
                graph.add_edge(predecessor, to_node, attr_dict=data)
            elif predecessor not in graph.pred[to_node] or key not in graph.edge[predecessor][to_node]:
                graph.add_edge(predecessor, to_node, key=key, **{RELATION: unqualified_edges[-1 - key]})

    graph.remove_node(from_node)


@pipeline.in_place_mutator
def collapse_nodes(graph, dict_of_sets_of_nodes):
    """Collapses all nodes in values to the key nodes, in place

    :param pybel.BELGraph graph: A BEL graph
    :param dict[tuple,set[tuple]] dict_of_sets_of_nodes: A dictionary of {node: set of nodes}
    """
    for key_node, value_nodes in dict_of_sets_of_nodes.items():
        for value_node in value_nodes:
            collapse_pair(graph, from_node=value_node, to_node=key_node)

    # Remove self edges
    for u, v, k in graph.edges(keys=True):
        if u == v:
            graph.remove_edge(u, v, k)


def build_central_dogma_collapse_dict(graph):
    """Builds a dictionary to direct the collapsing on the central dogma

    :param pybel.BELGraph graph: A BEL graph
    :return: A dictionary of {node: set of nodes}
    :rtype: dict[tuple,set[tuple]]
    """
    collapse_dict = defaultdict(set)
    r2p = {}

    for rna_node, protein_node, d in graph.edges_iter(data=True):
        if d[RELATION] != TRANSLATED_TO:
            continue

        collapse_dict[protein_node].add(rna_node)
        r2p[rna_node] = protein_node

    for gene_node, rna_node, d in graph.edges_iter(data=True):
        if d[RELATION] != TRANSCRIBED_TO:
            continue

        if rna_node in r2p:
            collapse_dict[r2p[rna_node]].add(gene_node)
        else:
            collapse_dict[rna_node].add(gene_node)

    return collapse_dict


def build_central_dogma_collapse_gene_dict(graph):
    """Builds a dictionary to direct the collapsing on the central dogma

    :param pybel.BELGraph graph: A BEL graph
    :return: A dictionary of {node: set of PyBEL node tuples}
    :rtype: dict[tuple,set[tuple]]
    """
    collapse_dict = defaultdict(set)
    r2g = {}

    for gene_node, rna_node, d in graph.edges_iter(data=True):
        if d[RELATION] != TRANSCRIBED_TO:
            continue

        collapse_dict[gene_node].add(rna_node)
        r2g[rna_node] = gene_node

    for rna_node, protein_node, d in graph.edges_iter(data=True):
        if d[RELATION] != TRANSLATED_TO:
            continue

        if rna_node not in r2g:
            raise ValueError('Should complete origin before running this function')

        collapse_dict[r2g[rna_node]].add(protein_node)

    return collapse_dict


@pipeline.in_place_mutator
def collapse_by_central_dogma(graph):
    """Collapses all nodes from the central dogma (GENE, RNA, PROTEIN) to PROTEIN, or most downstream possible entity,
    in place. This function wraps :func:`collapse_nodes` and :func:`build_central_dogma_collapse_dict`.

    :param pybel.BELGraph graph: A BEL graph
    
    Equivalent to:
    
    >>> collapse_nodes(graph, build_central_dogma_collapse_dict(graph)) 
    """
    collapse_dict = build_central_dogma_collapse_dict(graph)
    collapse_nodes(graph, collapse_dict)


@pipeline.in_place_mutator
def collapse_by_central_dogma_to_genes(graph):
    """Collapses all nodes from the central dogma (:data:`pybel.constants.GENE`, :data:`pybel.constants.RNA`, 
    :data:`pybel.constants.MIRNA`, and :data:`pybel.constants.PROTEIN`) to :data:`pybel.constants.GENE`, in-place. This 
    function wraps :func:`collapse_nodes` and :func:`build_central_dogma_collapse_gene_dict`.

    :param pybel.BELGraph graph: A BEL graph
    
    Equivalent to:
    
    >>> infer_central_dogma(graph)
    >>> collapse_nodes(graph, build_central_dogma_collapse_gene_dict(graph))
    """
    infer_central_dogma(graph)
    collapse_dict = build_central_dogma_collapse_gene_dict(graph)
    collapse_nodes(graph, collapse_dict)


@pipeline.mutator
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


@pipeline.in_place_mutator
def _collapse_variants_by_function(graph, func):
    """Collapses all of the given functions' variants' edges to their parents, in-place

    :param pybel.BELGraph graph: A BEL graph
    :param str func: A BEL function
    """
    for parent_node, variant_node, data in graph.edges(data=True):
        if data[RELATION] == HAS_VARIANT and graph.node[parent_node][FUNCTION] == func:
            collapse_pair(graph, from_node=variant_node, to_node=parent_node)


@pipeline.in_place_mutator
def collapse_protein_variants(graph):
    """Collapses all protein's variants' edges to their parents, in-place

    :param pybel.BELGraph graph: A BEL graph
    """
    _collapse_variants_by_function(graph, PROTEIN)


@pipeline.in_place_mutator
def collapse_gene_variants(graph):
    """Collapses all gene's variants' edges to their parents, in-place

    :param pybel.BELGraph graph: A BEL graph
    """
    _collapse_variants_by_function(graph, GENE)


@pipeline.in_place_mutator
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


@pipeline.in_place_mutator
def collapse_all_variants(graph):
    """Collapses all ``hasVariant`` edges to the parent node, in place
    
    :param pybel.BELGraph graph: A BEL Graph
    """
    for parent_node, variant_node, d in graph.edges(data=True):
        if d[RELATION] == HAS_VARIANT:
            collapse_pair(graph, from_node=variant_node, to_node=parent_node)


@pipeline.mutator
def collapse_all_variants_out_place(graph):
    """Collapses all ``hasVariant`` edges to the parent node, not in place

    :param pybel.BELGraph graph: A BEL Graph
    :rtype: pybel.BELGraph
    """
    result = graph.copy()
    _collapse_variants_by_function(result)
    return result


@pipeline.in_place_mutator
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


@pipeline.in_place_mutator
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


@pipeline.in_place_mutator
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


@pipeline.in_place_mutator
def collapse_orthologies_by_namespace(graph, from_namespace, to_namespace):
    """Collapses pairs of nodes with the given namespaces that have orthology relationships

    :param pybel.BELGraph graph: A BEL Graph
    :param str or iter[str] from_namespace: The namespace of the node to collapse
    :param str or iter[str] to_namespace: The namespace of the node to keep

    To collapse all MGI nodes to their HGNC orthologs, use:

    >>> collapse_orthologies_by_namespace('MGI', 'HGNC')
    """
    _collapse_edge_by_namespace(graph, from_namespace, to_namespace, ORTHOLOGOUS)


@pipeline.in_place_mutator
def collapse_entrez_to_hgnc(graph):
    """Collapses Entrez orthologies to HGNC.

    Implemented with :func:`collapse_orthologies_by_namespace`.

    :param pybel.BELGraph graph: A BEL graph
    """
    collapse_orthologies_by_namespace(graph, ['EGID', 'EG', 'ENTREZ'], 'HGNC')


@pipeline.in_place_mutator
def collapse_mgi_to_hgnc(graph):
    """Collapses MGI orthologies to HGNC

    Implemented with :func:`collapse_orthologies_by_namespace`.

    :param pybel.BELGraph graph: A BEL graph
    """
    collapse_orthologies_by_namespace(graph, ['MGI', 'MGIID'], 'HGNC')


@pipeline.in_place_mutator
def collapse_rgd_to_hgnc(graph):
    """Collapses RGD orthologies to HGNC

    Implemented with :func:`collapse_orthologies_by_namespace`.

    :param pybel.BELGraph graph: A BEL graph
    """
    collapse_orthologies_by_namespace(graph, ['RDG', 'RGDID'], 'HGNC')


@pipeline.in_place_mutator
def collapse_flybase_to_hgnc(graph):
    """Collapses FlyBase orthologies to HGNC

    Implemented with :func:`collapse_orthologies_by_namespace`.

    :param pybel.BELGraph graph: A BEL graph
    """
    collapse_orthologies_by_namespace(graph, 'FLYBASE', 'HGNC')


@pipeline.in_place_mutator
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


@pipeline.in_place_mutator
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


@pipeline.mutator
def collapse_to_protein_interactions(graph):
    """Collapses to a graph made of only causal gene/protein edges

    :param pybel.BELGraph graph: A BEL Graph
    """
    collapse_by_central_dogma_to_genes(graph)

    remove_filtered_edges(graph, invert_edge_filter(has_polarity))

    filtered_graph = graph.subgraph(get_nodes(graph, node_predicates=function_inclusion_filter_builder(GENE)))

    return filtered_graph
