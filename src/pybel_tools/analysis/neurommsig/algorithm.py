# -*- coding: utf-8 -*-

"""An implementation of the mechanism enrichment algorithm from Domingo-Fernández *et al.*, 2017."""

import itertools as itt
from collections import Counter

from pybel.constants import GENE
from ...filters.node_selection import get_nodes_by_function
from ...grouping import get_subgraphs_by_annotation
from ...mutation.inference import infer_central_dogma
from ...utils import calculate_betweenness_centality

__all__ = [
    'get_neurommsig_scores_stratified',
    'get_neurommsig_score',
]


def get_neurommsig_scores_stratified(graph, genes, annotation='Subgraph', ora_weight=None, hub_weight=None, top=None,
                                     topology_weight=None):
    """Runs the NeuroMMSig algorithm, stratified by a given annotation

    :param pybel.BELGraph graph: A BEL graph
    :param list[tuple] genes: A list of gene nodes
    :param str annotation: The annotation to use to stratify the graph to subgraphs
    :param Optional[float] ora_weight: The relative weight of the over-enrichment analysis score from
     :py:func:`neurommsig_gene_ora`. Defaults to 1.0.
    :param Optional[float] hub_weight: The relative weight of the hub analysis score from :py:func:`neurommsig_hubs`.
     Defaults to 1.0.
    :param Optional[float] top: The percentage of top genes to use as hubs. Defaults to 5% (0.05).
    :param Optional[float] topology_weight: The relative weight of the topolgical analysis core from
     :py:func:`neurommsig_topology`. Defaults to 1.0.
    :return: A dictionary from {annotation value: NeuroMMSig composite score}
    :rtype: dict[str, float]
    """
    return {
        annotation_value: get_neurommsig_score(subgraph, genes, ora_weight=ora_weight, hub_weight=hub_weight, top=top,
                                               topology_weight=topology_weight)
        for annotation_value, subgraph in get_subgraphs_by_annotation(graph, annotation=annotation).items()
    }


def get_neurommsig_score(graph, target_genes, ora_weight=None, hub_weight=None, top=None, topology_weight=None):
    """Calculates the composite NeuroMMSig Score for a given list of genes.

    :param pybel.BELGraph graph: A BEL graph
    :param list[tuple] target_genes: A list of gene nodes
    :param Optional[float] ora_weight: The relative weight of the over-enrichment analysis score from
     :py:func:`neurommsig_gene_ora`. Defaults to 1.0.
    :param Optional[float] hub_weight: The relative weight of the hub analysis score from :py:func:`neurommsig_hubs`.
     Defaults to 1.0.
    :param Optional[float] top: The percentage of top genes to use as hubs. Defaults to 5% (0.05).
    :param Optional[float] topology_weight: The relative weight of the topolgical analysis core from
     :py:func:`neurommsig_topology`. Defaults to 1.0.
    :return: The NeuroMMSig composite score
    :rtype: float
    """
    ora_weight = ora_weight or 1.0
    hub_weight = hub_weight or 1.0
    topology_weight = topology_weight or 1.0

    target_genes = list(target_genes)

    ora_score = neurommsig_gene_ora(graph, target_genes)
    hub_score = neurommsig_hubs(graph, target_genes, top=top)
    topology_score = neurommsig_topology(graph, target_genes)

    weighted_sum = ora_weight * ora_score + hub_weight * hub_score + topology_weight * topology_score
    total_weight = ora_weight + hub_weight + topology_weight
    return weighted_sum / total_weight


def neurommsig_gene_ora(graph, target_genes):
    """Calculates the percentage of target genes mappable to the graph
    
    Assume: graph central dogma inferred, collapsed to genes, collapsed variants 
    
    :param pybel.BELGraph graph: A BEL graph
    :param iter target_genes: An iterable of nodes
    :rtype: float
    """
    infer_central_dogma(graph)
    graph_genes = set(get_nodes_by_function(graph, GENE))
    return len(graph_genes.intersection(target_genes)) / len(graph_genes)


def neurommsig_hubs(graph, target_genes, top=None):
    """Calculates the percentage of target genes mappable to the graph
    
    Assume: graph central dogma inferred, collapsed to genes, collapsed variants, graph has more than 20 nodes
    
    :param pybel.BELGraph graph: A BEL graph
    :param iter[tuple] target_genes: A list of nodes
    :param Optional[float] top: The percentage of top genes to use as hubs. Defaults to 5% (0.05).
    :rtype: float
    """
    top = top or 0.05

    if graph.number_of_nodes() < 20:
        raise ValueError('Graph has less than 20 nodes')

    graph_genes = set(get_nodes_by_function(graph, GENE))

    bc = Counter({
        k: v
        for k, v in calculate_betweenness_centality(graph).items()
        if k in graph_genes
    })

    # TODO consider continuous analog with weighting by percentile
    n = int(0.5 + len(graph_genes) * top)
    unnormalized_sum = sum(
        node in target_genes
        for node in bc.most_common(n)
    )

    return unnormalized_sum / n


def neurommsig_topology(graph, nodes):
    """Calculates the node neighbor score for a given list of nodes.
    
    -  Doesn't consider self loops
    
    :param pybel.BELGraph graph: A BEL graph
    :param list[tuple] nodes: A list of nodes
    :rtype: float
    
    .. math::
        
         \frac{\sum_i^n N_G[i]}{n*(n-1)}
    """
    nodes = list(nodes)
    n = len(nodes)

    unnormalized_sum = sum(
        u in graph[v]
        for u, v in itt.product(nodes, repeat=2)
        if u != v
    )

    return unnormalized_sum / (n * (n - 1))
