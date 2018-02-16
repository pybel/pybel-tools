# -*- coding: utf-8 -*-

"""This module contains convenient functions for removing nodes/edges that are returned from selection functions"""

from pybel.constants import BIOPROCESS, PATHOLOGY
from pybel.struct.filters import get_nodes
from pybel.struct.filters.edge_filters import filter_edges
from pybel.struct.filters.edge_predicates import is_associative_relation
from .. import pipeline
from ..filters.node_selection import function_inclusion_filter_builder
from ..selection.leaves import get_gene_leaves, get_rna_leaves
from ..selection.utils import get_leaves_by_type
from ..summary.edge_summary import get_inconsistent_edges
from ..utils import all_edges_iter

__all__ = [
    'remove_filtered_edges',
    'remove_leaves_by_type',
    'prune_central_dogma',
    'remove_inconsistent_edges',
    'remove_pathologies',
    'remove_biological_processes',
    'remove_associations',
]


@pipeline.in_place_mutator
def remove_filtered_edges(graph, edge_filters):
    """Removes edges passing the given edge filters

    :param pybel.BELGraph graph: A BEL graph
    :param edge_filters: An edge filter or list of edge filters (graph, node, node, key, data)-> bool
    :type edge_filters: types.FunctionType or iter[types.FunctionType]
    :return: 
    """
    edges = list(filter_edges(graph, edge_filters))
    graph.remove_edges_from(edges)


@pipeline.in_place_mutator
def remove_leaves_by_type(graph, function=None, prune_threshold=1):
    """Removes all nodes in graph (in-place) with only a connection to one node. Useful for gene and RNA.
    Allows for optional filter by function type.


    :param pybel.BELGraph graph: A BEL graph
    :param function: If set, filters by the node's function from :mod:`pybel.constants` like
                    :data:`pybel.constants.GENE`, :data:`pybel.constants.RNA`, :data:`pybel.constants.PROTEIN`, or
                    :data:`pybel.constants.BIOPROCESS`
    :type function: str
    :param prune_threshold: Removes nodes with less than or equal to this number of connections. Defaults to :code:`1`
    :type prune_threshold: int
    """
    nodes = list(get_leaves_by_type(graph, function=function, prune_threshold=prune_threshold))
    graph.remove_nodes_from(nodes)


@pipeline.in_place_mutator
def prune_central_dogma(graph):
    """Prunes genes, then RNA, in place

    :param pybel.BELGraph graph: A BEL graph
    """
    gene_leaves = list(get_gene_leaves(graph))
    graph.remove_nodes_from(gene_leaves)

    rna_leaves = list(get_rna_leaves(graph))
    graph.remove_nodes_from(rna_leaves)


@pipeline.in_place_mutator
def remove_inconsistent_edges(graph):
    """Remove all edges between node paris with consistent edges.

    This is the all-or-nothing approach. It would be better to do more careful investigation of the evidences during
    curation.

    :param pybel.BELGraph graph: A BEL graph
    """
    for u, v in get_inconsistent_edges(graph):
        edges = list(all_edges_iter(graph, u, v))
        graph.remove_edges_from(edges)


@pipeline.in_place_mutator
def remove_pathologies(graph):
    """Remove pathology nodes

    :param pybel.BELGraph graph: A BEL graph
    """
    nodes = get_nodes(graph, function_inclusion_filter_builder(PATHOLOGY))
    graph.remove_nodes_from(nodes)


@pipeline.in_place_mutator
def remove_biological_processes(graph):
    """Remove biological process nodes

    :param pybel.BELGraph graph: A BEL graph
    """
    nodes = get_nodes(graph, function_inclusion_filter_builder(BIOPROCESS))
    graph.remove_nodes_from(nodes)


@pipeline.in_place_mutator
def remove_associations(graph):
    """Removes all associative relationships from the graph

    :param pybel.BELGraph graph: A BEL graph
    """
    remove_filtered_edges(graph, is_associative_relation)
