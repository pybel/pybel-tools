# -*- coding: utf-8 -*-

from collections import defaultdict

from pybel.constants import *
from pybel.struct.filters.node_filters import concatenate_node_filters
from ..utils import check_has_annotation

__all__ = [
    'group_nodes_by_annotation',
    'average_node_annotation',
    'group_nodes_by_annotation_filtered'
]


def group_nodes_by_annotation(graph, annotation='Subgraph'):
    """Groups the nodes occurring in edges by the given annotation

    :param pybel.BELGraph graph: A BEL graph
    :param annotation: An annotation to use to group edges
    :type annotation: str
    :return: dict of sets of BELGraph nodes
    :rtype: dict
    """
    result = defaultdict(set)

    for u, v, d in graph.edges_iter(data=True):
        if not check_has_annotation(d, annotation):
            continue

        result[d[ANNOTATIONS][annotation]].add(u)
        result[d[ANNOTATIONS][annotation]].add(v)

    return dict(result)


def average_node_annotation(graph, key, annotation='Subgraph', aggregator=None):
    """Groups graph into subgraphs and assigns each subgraph a score based on the average of all nodes values
    for the given node key

    :param pybel.BELGraph graph: A BEL graph
    :param key: The key in the node data dictionary representing the experimental data
    :type key: str
    :param annotation: A BEL annotation to use to group nodes
    :type annotation: str
    :param aggregator: A function from list of values -> aggregate value. Defaults to taking the average of a list of
                       floats.
    :type aggregator: lambda
    """

    if aggregator is None:
        def aggregator(x):
            """Calculates the average"""
            return sum(x) / len(x)

    result = {}
    grouping = group_nodes_by_annotation(graph, annotation)
    for subgraph, nodes in grouping.items():
        values = [graph.node[node][key] for node in nodes if key in graph.node[node]]
        result[subgraph] = aggregator(values)
    return result


def group_nodes_by_annotation_filtered(graph, node_filters=None, annotation='Subgraph'):
    """Groups the nodes occurring in edges by the given annotation, with a node filter applied

    :param pybel.BELGraph graph: A BEL graph
    :param node_filters: A predicate or list of predicates (graph, node) -> bool
    :type node_filters: types.FunctionType or iter[types.FunctionType]
    :param annotation: The annotation to use for grouping
    :return: A dictionary of {annotation value: set of nodes}
    :rtype: dict[str,set[tuple]]
    """
    node_filter = concatenate_node_filters(node_filters)
    return {
        key: {
            node
            for node in nodes
            if node_filter(graph, node)
        }
        for key, nodes in group_nodes_by_annotation(graph, annotation).items()
    }


def get_mapped(graph, namespace, names):
    """ Returns defaultdict with keys: nodes that match the namespace and in names and values other nodes (complexes, variants, orthologous...) or this node.
    
    :param pybel.BELGraph graph: A BEL graph
    :param str namespace: Namespace
    :param list[str] or set[str] names: List or set of values from which we want to map nodes from
    :return:
    """
    mapped_nodes = defaultdict(set)

    for u, v, d in graph.edges_iter(data=True):

        if d[RELATION] in {HAS_MEMBER, HAS_COMPONENT} and NAMESPACE in graph.node[v] and graph.node[v][
            NAMESPACE] == namespace and graph.node[v][NAME] in names:
            mapped_nodes[v].add(u)

        elif d[RELATION] == HAS_VARIANT and NAMESPACE in graph.node[u] and graph.node[u][NAMESPACE] == namespace and \
                        graph.node[u][NAME] in names:
            mapped_nodes[u].add(v)

        elif d[RELATION] == ORTHOLOGOUS and NAMESPACE in graph.node[u] and graph.node[u][NAMESPACE] == namespace and \
                        graph.node[u][NAME] in names:
            mapped_nodes[u].add(v)

    mapped_nodes.default_factory = None

    return mapped_nodes
