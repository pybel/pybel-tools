# -*- coding: utf-8 -*-

"""Collapse functions to supplement :mod:`pybel.struct.mutation.collapse`."""

import itertools as itt
import logging
from collections import defaultdict

import networkx as nx
from pybel import BELGraph
from pybel.constants import EQUIVALENT_TO, GENE, HAS_VARIANT, NAME, NAMESPACE, ORTHOLOGOUS, PROTEIN, RELATION
from pybel.dsl import BaseEntity, Gene, Protein
from pybel.struct.filters import build_relation_predicate, filter_edges, has_polarity
from pybel.struct.filters.typing import EdgePredicates
from pybel.struct.mutation import collapse_nodes, collapse_pair, collapse_to_genes, get_subgraph_by_edge_filter
from pybel.struct.pipeline import in_place_transformation, transformation
from pybel.typing import Strings
from tqdm import tqdm

from ..filters.edge_filters import build_source_namespace_filter, build_target_namespace_filter
from ..summary.edge_summary import pair_is_consistent

__all__ = [
    'collapse_nodes',
    'rewire_variants_to_genes',
    'collapse_gene_variants',
    'collapse_protein_variants',
    'collapse_consistent_edges',
    'collapse_equivalencies_by_namespace',
    'collapse_orthologies_by_namespace',
    'collapse_to_protein_interactions',
    'collapse_nodes_with_same_names',
]

log = logging.getLogger(__name__)


@in_place_transformation
def collapse_protein_variants(graph: BELGraph) -> None:
    """Collapse all protein's variants' edges to their parents, in-place."""
    _collapse_variants_by_function(graph, PROTEIN)


@in_place_transformation
def collapse_gene_variants(graph: BELGraph) -> None:
    """Collapse all gene's variants' edges to their parents, in-place."""
    _collapse_variants_by_function(graph, GENE)


def _collapse_variants_by_function(graph: BELGraph, func: str) -> None:
    """Collapse all of the given functions' variants' edges to their parents, in-place."""
    for parent_node, variant_node, data in graph.edges(data=True):
        if data[RELATION] == HAS_VARIANT and parent_node.function == func:
            collapse_pair(graph, from_node=variant_node, to_node=parent_node)


@in_place_transformation
def rewire_variants_to_genes(graph: BELGraph) -> None:
    """Find all protein variants that are pointing to a gene and not a protein and fixes them by changing their
    function to be :data:`pybel.constants.GENE`, in place

    A use case is after running :func:`collapse_to_genes`.
    """
    mapping = {}

    for node in graph:
        if not isinstance(node, Protein) or not node.variants:
            continue

        mapping[node] = Gene(
            name=node.name,
            namespace=node.namespace,
            identifier=node.identifier,
            variants=node.variants,
        )

    nx.relabel_nodes(graph, mapping, copy=False)


def _collapse_edge_passing_predicates(graph: BELGraph, edge_predicates: EdgePredicates = None) -> None:
    """Collapse all edges passing the given edge predicates."""
    for u, v, _ in filter_edges(graph, edge_predicates=edge_predicates):
        collapse_pair(graph, survivor=u, victim=v)


def _collapse_edge_by_namespace(graph: BELGraph,
                                victim_namespaces: Strings,
                                survivor_namespaces: str,
                                relations: Strings) -> None:
    """Collapse pairs of nodes with the given namespaces that have the given relationship.

    :param graph: A BEL Graph
    :param victim_namespaces: The namespace(s) of the node to collapse
    :param survivor_namespaces: The namespace of the node to keep
    :param relations: The relation(s) to search
    """
    relation_filter = build_relation_predicate(relations)
    source_namespace_filter = build_source_namespace_filter(victim_namespaces)
    target_namespace_filter = build_target_namespace_filter(survivor_namespaces)

    edge_predicates = [
        relation_filter,
        source_namespace_filter,
        target_namespace_filter
    ]

    _collapse_edge_passing_predicates(graph, edge_predicates=edge_predicates)


@in_place_transformation
def collapse_equivalencies_by_namespace(graph: BELGraph, victim_namespace: Strings, survivor_namespace: str) -> None:
    """Collapse pairs of nodes with the given namespaces that have equivalence relationships.
    
    :param graph: A BEL graph
    :param victim_namespace: The namespace(s) of the node to collapse
    :param survivor_namespace: The namespace of the node to keep

    To convert all ChEBI names to InChI keys, assuming there are appropriate equivalence relations between nodes with
    those namespaces:
    
    >>> collapse_equivalencies_by_namespace(graph, 'CHEBI', 'CHEBIID')
    >>> collapse_equivalencies_by_namespace(graph, 'CHEBIID', 'INCHI')
    """
    _collapse_edge_by_namespace(graph, victim_namespace, survivor_namespace, EQUIVALENT_TO)


@in_place_transformation
def collapse_orthologies_by_namespace(graph: BELGraph, victim_namespace: Strings, survivor_namespace: str) -> None:
    """Collapse pairs of nodes with the given namespaces that have orthology relationships.

    :param graph: A BEL Graph
    :param victim_namespace: The namespace(s) of the node to collapse
    :param survivor_namespace: The namespace of the node to keep

    To collapse all MGI nodes to their HGNC orthologs, use:
    >>> collapse_orthologies_by_namespace('MGI', 'HGNC')


    To collapse collapse both MGI and RGD nodes to their HGNC orthologs, use:
    >>> collapse_orthologies_by_namespace(['MGI', 'RGD'], 'HGNC')
    """
    _collapse_edge_by_namespace(graph, victim_namespace, survivor_namespace, ORTHOLOGOUS)


@in_place_transformation
def collapse_entrez_to_hgnc(graph: BELGraph):
    """Collapse Entrez equivalences to HGNC."""
    collapse_equivalencies_by_namespace(graph, ['EGID', 'EG', 'ENTREZ'], 'HGNC')


@in_place_transformation
def collapse_mgi_to_hgnc(graph: BELGraph):
    """Collapse MGI orthologies to HGNC."""
    collapse_orthologies_by_namespace(graph, ['MGI', 'MGIID'], 'HGNC')


@in_place_transformation
def collapse_rgd_to_hgnc(graph: BELGraph):
    """Collapse RGD orthologies to HGNC."""
    collapse_orthologies_by_namespace(graph, ['RGD', 'RGDID'], 'HGNC')


@in_place_transformation
def collapse_flybase_to_hgnc(graph: BELGraph):
    """Collapse FlyBase orthologies to HGNC."""
    collapse_orthologies_by_namespace(graph, 'FLYBASE', 'HGNC')


@in_place_transformation
def collapse_entrez_equivalencies(graph: BELGraph):
    """Collapse all equivalence edges away from Entrez. Assumes well formed, 2-way equivalencies."""
    relation_filter = build_relation_predicate(EQUIVALENT_TO)
    source_namespace_filter = build_source_namespace_filter(['EGID', 'EG', 'ENTREZ'])

    edge_predicates = [
        relation_filter,
        source_namespace_filter,
    ]

    _collapse_edge_passing_predicates(graph, edge_predicates=edge_predicates)


@in_place_transformation
def collapse_consistent_edges(graph: BELGraph):
    """Collapse consistent edges together.

    .. warning:: This operation doesn't preserve evidences or other annotations
    """
    for u, v in graph.edges():
        relation = pair_is_consistent(graph, u, v)

        if not relation:
            continue

        edges = [(u, v, k) for k in graph[u][v]]
        graph.remove_edges_from(edges)
        graph.add_edge(u, v, attr_dict={RELATION: relation})


@transformation
def collapse_to_protein_interactions(graph: BELGraph) -> BELGraph:
    """Collapse to a graph made of only causal gene/protein edges."""
    rv: BELGraph = graph.copy()

    collapse_to_genes(rv)

    def is_edge_ppi(_: BELGraph, u: BaseEntity, v: BaseEntity, __: str) -> bool:
        """Check if an edge is a PPI."""
        return isinstance(u, Gene) and isinstance(v, Gene)

    return get_subgraph_by_edge_filter(rv, edge_predicates=[has_polarity, is_edge_ppi])


@in_place_transformation
def collapse_nodes_with_same_names(graph: BELGraph) -> None:
    """Collapse all nodes with the same name, merging namespaces by picking first alphabetical one."""
    survivor_mapping = defaultdict(set) # Collapse mapping dict
    victims = set() # Things already mapped while iterating

    it = tqdm(itt.combinations(graph, r=2), total=graph.number_of_nodes() * (graph.number_of_nodes() - 1) / 2)
    for a, b in it:
        if b in victims:
            continue

        a_name, b_name = a.get(NAME), b.get(NAME)
        if not a_name or not b_name or a_name.lower() != b_name.lower():
            continue

        if a.keys() != b.keys():  # not same version (might have variants)
            continue

        # Ensure that the values in the keys are also the same
        for k in set(a.keys()) - {NAME, NAMESPACE}:
            if a[k] != b[k]:  # something different
                continue

        survivor_mapping[a].add(b)
        # Keep track of things that has been already mapped
        victims.add(b)

    collapse_nodes(graph, survivor_mapping)
