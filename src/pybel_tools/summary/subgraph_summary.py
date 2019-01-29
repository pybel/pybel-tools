# -*- coding: utf-8 -*-

"""This module contains functions that handle and summarize subgraphs of graphs."""

import itertools as itt
from collections import defaultdict
from operator import itemgetter
from typing import Counter, Iterable, List, Mapping, Set, Tuple, Union

from pybel import BELGraph
from pybel.constants import ANNOTATIONS
from pybel.dsl import BaseEntity
from pybel.struct.filters.edge_predicates import edge_has_annotation
from pybel.struct.filters.typing import NodePredicate
from ..selection.group_nodes import group_nodes_by_annotation, group_nodes_by_annotation_filtered
from ..utils import calculate_tanimoto_set_distances, count_dict_values

__all__ = [
    'count_subgraph_sizes',
    'calculate_subgraph_edge_overlap',
    'summarize_subgraph_edge_overlap',
    'rank_subgraph_by_node_filter',
    'summarize_subgraph_node_overlap',
]


def count_subgraph_sizes(graph: BELGraph, annotation: str = 'Subgraph') -> Counter[int]:
    """Count the number of nodes in each subgraph induced by an annotation.

    :param annotation: The annotation to group by and compare. Defaults to 'Subgraph'
    :return: A dictionary from {annotation value: number of nodes}
    """
    return count_dict_values(group_nodes_by_annotation(graph, annotation))


EdgeSet = Set[Tuple[BaseEntity, BaseEntity]]


def calculate_subgraph_edge_overlap(
        graph: BELGraph,
        annotation: str = 'Subgraph'
) -> Tuple[
    Mapping[str, EdgeSet],
    Mapping[str, Mapping[str, EdgeSet]],
    Mapping[str, Mapping[str, EdgeSet]],
    Mapping[str, Mapping[str, float]],
]:
    """Build a DatafFame to show the overlap between different sub-graphs.

    Options:
    1. Total number of edges overlap (intersection)
    2. Percentage overlap (tanimoto similarity)

    :param graph: A BEL graph
    :param annotation: The annotation to group by and compare. Defaults to 'Subgraph'
    :return: {subgraph: set of edges}, {(subgraph 1, subgraph2): set of intersecting edges},
            {(subgraph 1, subgraph2): set of unioned edges}, {(subgraph 1, subgraph2): tanimoto similarity},
    """
    sg2edge = defaultdict(set)

    for u, v, d in graph.edges(data=True):
        if not edge_has_annotation(d, annotation):
            continue
        sg2edge[d[ANNOTATIONS][annotation]].add((u, v))

    subgraph_intersection = defaultdict(dict)
    subgraph_union = defaultdict(dict)
    result = defaultdict(dict)

    for sg1, sg2 in itt.product(sg2edge, repeat=2):
        subgraph_intersection[sg1][sg2] = sg2edge[sg1] & sg2edge[sg2]
        subgraph_union[sg1][sg2] = sg2edge[sg1] | sg2edge[sg2]
        result[sg1][sg2] = len(subgraph_intersection[sg1][sg2]) / len(subgraph_union[sg1][sg2])

    return sg2edge, subgraph_intersection, subgraph_union, result


def summarize_subgraph_edge_overlap(graph: BELGraph, annotation: str = 'Subgraph') -> Mapping[str, Mapping[str, float]]:
    """Return a similarity matrix between all subgraphs (or other given annotation).

    :param annotation: The annotation to group by and compare. Defaults to :code:`"Subgraph"`
    :return: A similarity matrix in a dict of dicts
    :rtype: dict
    """
    _, _, _, subgraph_overlap = calculate_subgraph_edge_overlap(graph, annotation)
    return subgraph_overlap


def summarize_subgraph_node_overlap(graph: BELGraph, node_predicates=None, annotation: str = 'Subgraph'):
    """Calculate the subgraph similarity tanimoto similarity in nodes passing the given filter.

    Provides an alternate view on subgraph similarity, from a more node-centric view
    """
    r1 = group_nodes_by_annotation_filtered(graph, node_predicates=node_predicates, annotation=annotation)
    return calculate_tanimoto_set_distances(r1)


def rank_subgraph_by_node_filter(graph: BELGraph,
                                 node_predicates: Union[NodePredicate, Iterable[NodePredicate]],
                                 annotation: str = 'Subgraph',
                                 reverse: bool = True,
                                 ) -> List[Tuple[str, int]]:
    """Rank sub-graphs by which have the most nodes matching an given filter.

    A use case for this function would be to identify which subgraphs contain the most differentially expressed
    genes.

    >>> from pybel import from_pickle
    >>> from pybel.constants import GENE
    >>> from pybel_tools.integration import overlay_type_data
    >>> from pybel_tools.summary import rank_subgraph_by_node_filter
    >>> import pandas as pd
    >>> graph = from_pickle('~/dev/bms/aetionomy/alzheimers.gpickle')
    >>> df = pd.read_csv('~/dev/bananas/data/alzheimers_dgxp.csv', columns=['Gene', 'log2fc'])
    >>> data = {gene: log2fc for _, gene, log2fc in df.itertuples()}
    >>> overlay_type_data(graph, data, 'log2fc', GENE, 'HGNC', impute=0.0)
    >>> results = rank_subgraph_by_node_filter(graph, lambda g, n: 1.3 < abs(g[n]['log2fc']))
    """
    r1 = group_nodes_by_annotation_filtered(graph, node_predicates=node_predicates, annotation=annotation)
    r2 = count_dict_values(r1)
    # TODO use instead: r2.most_common()
    return sorted(r2.items(), key=itemgetter(1), reverse=reverse)
