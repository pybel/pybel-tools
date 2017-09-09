# -*- coding: utf-8 -*-

import logging

from collections import defaultdict

from pybel.constants import *
from .deletion import prune_central_dogma
from .inference import infer_central_dogma
from .. import pipeline
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
    'opening_on_central_dogma',
    'collapse_consistent_edges',
]

log = logging.getLogger(__name__)


@pipeline.in_place_mutator
def collapse_pair(graph, survivor, synonym):
    """Rewires all edges from the synonymous node to the survivor node, then deletes the synonymous node.
    
    Does not keep edges between the two nodes.
    
    :param pybel.BELGraph graph: A BEL graph
    :param tuple survivor: The BEL node to collapse all edges on the synonym to
    :param tuple synonym: The BEL node to collapse into the surviving node
    """
    for successor in graph.successors_iter(synonym):
        if successor == survivor:
            continue

        for key, data in graph.edge[synonym][successor].items():
            if key >= 0:
                graph.add_edge(survivor, successor, attr_dict=data)
            elif successor not in graph.edge[survivor] or key not in graph.edge[survivor][successor]:
                graph.add_edge(survivor, successor, key=key, **{RELATION: unqualified_edges[-1 - key]})

    for predecessor in graph.predecessors_iter(synonym):
        if predecessor == survivor:
            continue

        for key, data in graph.edge[predecessor][synonym].items():
            if key >= 0:
                graph.add_edge(predecessor, survivor, attr_dict=data)
            elif predecessor not in graph.pred[survivor] or key not in graph.edge[predecessor][survivor]:
                graph.add_edge(predecessor, survivor, key=key, **{RELATION: unqualified_edges[-1 - key]})

    graph.remove_node(synonym)


@pipeline.in_place_mutator
def collapse_nodes(graph, dict_of_sets_of_nodes):
    """Collapses all nodes in values to the key nodes, in place

    :param pybel.BELGraph graph: A BEL graph
    :param dict[tuple,set[tuple]] dict_of_sets_of_nodes: A dictionary of {node: set of nodes}
    """
    log.debug('collapsing %d groups', len(dict_of_sets_of_nodes))

    for key_node, value_nodes in dict_of_sets_of_nodes.items():
        for value_node in value_nodes:
            collapse_pair(graph, key_node, value_node)

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
def _collapse_variants_by_function(graph, function):
    """Collapses all of the given functions' variants' edges to their parents, in-place

    :param pybel.BELGraph graph: A BEL graph
    :param str function: A BEL function
    """
    for parent_node, variant_node, d in graph.edges(data=True):
        if d[RELATION] == HAS_VARIANT and graph.node[parent_node][FUNCTION] == function:
            collapse_pair(graph, parent_node, variant_node)


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
    for u, v, d in graph.edges(data=True):
        if d[RELATION] == HAS_VARIANT:
            collapse_pair(graph, u, v)


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
def opening_on_central_dogma(graph):
    """Infers central dogmatic relations with :func:`infer_central_dogma` then successively prunes gene leaves then
    RNA leaves with :func:`prune_central_dogma` to connect disparate elements in a knowledge assembly

    :param pybel.BELGraph graph: A BEL graph

    Equivalent to:

    >>> infer_central_dogma(graph)
    >>> prune_central_dogma(graph)

    """
    infer_central_dogma(graph)
    prune_central_dogma(graph)


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


@pipeline.in_place_mutator
def collapse_namespace(graph, from_namespace, to_namespace):
    """Collapses pairs of nodes that have equivalence relationships to the target namespace
    
    :param pybel.BELGraph graph: A BEL graph
    :param str from_namespace: 
    :param str to_namespace: 
    
    
    To convert all ChEBI names to InChI keys:
    
    >>> collapse_namespace(graph, 'CHEBI', 'CHEBIID')
    >>> collapse_namespace(graph, 'CHEBIID', 'INCHI')
    """
    for u, v, d in graph.edges(data=True):
        if d[RELATION] != EQUIVALENT_TO:
            continue

        if NAMESPACE not in graph.node[u] or graph.node[u][NAMESPACE] != from_namespace:
            continue

        if NAMESPACE not in graph.node[v] or graph.node[v][NAMESPACE] != to_namespace:
            continue

        collapse_pair(graph, u, v)


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
