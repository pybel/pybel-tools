# -*- coding: utf-8 -*-

"""An implementation of the mechanism enrichment algorithm from Domingo-Fernández *et al.*, 2017."""

import itertools as itt
import logging
from collections import Counter

from pybel.constants import GENE

from ...filters.node_selection import get_nodes_by_function
from ...grouping import get_subgraphs_by_annotation
from ...mutation import collapse_all_variants, collapse_by_central_dogma_to_genes, infer_central_dogma
from ...pipeline import Pipeline
from ...utils import calculate_betweenness_centality

__all__ = [
    'neurommsig_graph_preprocessor',
    'get_neurommsig_scores_prestratified',
    'get_neurommsig_scores',
    'get_neurommsig_score',
]

log = logging.getLogger(__name__)

neurommsig_graph_preprocessor = Pipeline.from_functions([
    infer_central_dogma,
    collapse_by_central_dogma_to_genes,
    collapse_all_variants,
])


def get_neurommsig_scores_prestratified(it, genes, ora_weight=None, hub_weight=None, top=None, topology_weight=None):
    """Takes a graph stratification and runs neurommsig on each

    :param iter[tuple[str,pybel.BELGraph]] it: A pre-stratified set of graphs
    :param list[tuple] genes: A list of gene nodes
    :param Optional[float] ora_weight: The relative weight of the over-enrichment analysis score from
     :py:func:`neurommsig_gene_ora`. Defaults to 1.0.
    :param Optional[float] hub_weight: The relative weight of the hub analysis score from :py:func:`neurommsig_hubs`.
     Defaults to 1.0.
    :param Optional[float] top: The percentage of top genes to use as hubs. Defaults to 5% (0.05).
    :param Optional[float] topology_weight: The relative weight of the topolgical analysis core from
     :py:func:`neurommsig_topology`. Defaults to 1.0.
    :param bool preprocess: If true, preprocess the graph.
    :return: A dictionary from {annotation value: NeuroMMSig composite score}
    :rtype: Optional[dict[str, float]]

    Pre-processing steps:

    1. Infer the central dogma with :func:``
    2. Collapse all proteins, RNAs and miRNAs to genes with :func:``
    3. Collapse variants to genes with :func:``
    """
    rv = {}

    for annotation_value, subgraph in it:
        score = get_neurommsig_score(subgraph, genes, ora_weight=ora_weight, hub_weight=hub_weight, top=top,
                                     topology_weight=topology_weight)
        rv[annotation_value] = score

    return rv


def get_neurommsig_scores(graph, genes, annotation='Subgraph', ora_weight=None, hub_weight=None, top=None,
                          topology_weight=None, preprocess=False):
    """Preprocesses the graph, stratifies by the given annotation, then runs the NeuroMMSig algorithm on each.

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
    :param bool preprocess: If true, preprocess the graph.
    :return: A dictionary from {annotation value: NeuroMMSig composite score}
    :rtype: Optional[dict[str, float]]

    Pre-processing steps:

    1. Infer the central dogma with :func:``
    2. Collapse all proteins, RNAs and miRNAs to genes with :func:``
    3. Collapse variants to genes with :func:``
    """
    if preprocess:
        graph = neurommsig_graph_preprocessor.run(graph)

    if not any(gene in graph for gene in genes):
        log.debug('no genes mapping to graph')
        return

    it = get_subgraphs_by_annotation(graph, annotation=annotation).items()

    return get_neurommsig_scores_prestratified(it, genes, ora_weight=ora_weight, hub_weight=hub_weight, top=top,
                                               topology_weight=topology_weight)


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
        log.debug('Graph has less than 20 nodes')
        return 0.0

    graph_genes = set(get_nodes_by_function(graph, GENE))

    bc = Counter({
        node: betweenness_centrality
        for node, betweenness_centrality in calculate_betweenness_centality(graph).items()
        if node in graph_genes
    })

    # TODO consider continuous analog with weighting by percentile
    n = int(len(graph_genes) * top)

    if n < 1:
        n = 1

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

    if n <= 1:
        # log.debug('')
        return 0.0

    unnormalized_sum = sum(
        u in graph[v]
        for u, v in itt.product(nodes, repeat=2)
        if v in graph and u != v
    )

    return unnormalized_sum / (n * (n - 1.0))
