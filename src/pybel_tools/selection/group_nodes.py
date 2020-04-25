# -*- coding: utf-8 -*-

"""Node grouping utilities."""

from collections import defaultdict
from typing import Callable, Iterable, List, Mapping, Optional, Set, TypeVar

from pybel import BELGraph, BaseEntity
from pybel.constants import ANNOTATIONS, HAS_VARIANT, IS_A, ORTHOLOGOUS, PART_OF, RELATION
from pybel.dsl import BaseConcept
from pybel.struct.filters import concatenate_node_predicates
from pybel.struct.filters.edge_predicates import edge_has_annotation
from pybel.struct.filters.typing import NodePredicates
from ..utils import group_as_sets

__all__ = [
    'group_nodes_by_annotation',
    'average_node_annotation',
    'group_nodes_by_annotation_filtered',
    'get_mapped_nodes',
]

X = TypeVar('X')


def group_nodes_by_annotation(
    graph: BELGraph,
    annotation: str = 'Subgraph',
) -> Mapping[str, Set[BaseEntity]]:
    """Group the nodes occurring in edges by the given annotation."""
    result = defaultdict(set)

    for u, v, d in graph.edges(data=True):
        if not edge_has_annotation(d, annotation):
            continue

        result[d[ANNOTATIONS][annotation]].add(u)
        result[d[ANNOTATIONS][annotation]].add(v)

    return dict(result)


def average_node_annotation(
    graph: BELGraph,
    key: str,
    annotation: str = 'Subgraph',
    aggregator: Optional[Callable[[Iterable[X]], X]] = None,
) -> Mapping[str, X]:
    """Group a graph into sub-graphs and calculate an aggregate score for all nodes in each.

    :param graph: A BEL graph
    :param key: The key in the node data dictionary representing the experimental data
    :param annotation: A BEL annotation to use to group nodes
    :param aggregator: A function from list of values -> aggregate value. Defaults to taking the average of a list of
                       floats.
    """
    if aggregator is None:
        def aggregator(x: List[float]) -> float:
            """Calculate the average."""
            return sum(x) / len(x)

    grouped_nodes = group_nodes_by_annotation(graph, annotation)

    return {
        subgraph: aggregator([
            graph.nodes[node][key]
            for node in nodes
            if key in graph.nodes[node]
        ])
        for subgraph, nodes in grouped_nodes.items()
    }


def group_nodes_by_annotation_filtered(
    graph: BELGraph,
    node_predicates: NodePredicates,
    annotation: str = 'Subgraph',
) -> Mapping[str, Set[BaseEntity]]:
    """Group the nodes occurring in edges by the given annotation, with a node filter applied.

    :param graph: A BEL graph
    :param node_predicates: A predicate or list of predicates (graph, node) -> bool
    :param annotation: The annotation to use for grouping
    :return: A dictionary of {annotation value: set of nodes}
    """
    if node_predicates is None:
        raise ValueError('Just use group_nodes_by_annotation() if you do not have a node filter')

    node_predicate = concatenate_node_predicates(node_predicates)

    return {
        key: {
            node
            for node in nodes
            if node_predicate(graph, node)
        }
        for key, nodes in group_nodes_by_annotation(graph, annotation).items()
    }


def get_mapped_nodes(
    graph: BELGraph,
    namespace: str,
    names: Iterable[str],
) -> Mapping[BaseEntity, Set[BaseEntity]]:
    """Get the nodes mapped to this node's concept.

    Returns a dict with keys: nodes that match the namespace and in
    names and values other nodes (complexes, variants, orthologous...)
    or this node.

    :param graph: A BEL graph
    :param namespace: The namespace to search
    :param names: List or set of values from which we want to map nodes from
    :return: Main node to variants/groups.
    """
    names = {n.lower() for n in names}
    namespace = namespace.lower()

    return group_as_sets(
        (parent_node, mapped_node)
        for parent_node, mapped_node, d in graph.edges(data=True)
        if (
            isinstance(parent_node, BaseConcept)
            and parent_node.namespace.lower() == namespace
            and parent_node.name.lower() in names
            and d[RELATION] in {IS_A, PART_OF, HAS_VARIANT, ORTHOLOGOUS}
        )
    )
