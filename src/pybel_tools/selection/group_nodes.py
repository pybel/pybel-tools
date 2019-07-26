# -*- coding: utf-8 -*-

"""Node grouping utilities."""

from collections import defaultdict
from typing import Callable, Iterable, List, Mapping, Optional, Set, TypeVar

from pybel import BELGraph, BaseEntity
from pybel.constants import ANNOTATIONS, HAS_COMPONENT, HAS_MEMBER, HAS_VARIANT, NAME, NAMESPACE, ORTHOLOGOUS, RELATION
from pybel.struct.filters.edge_predicates import edge_has_annotation
from pybel.struct.filters.node_filters import concatenate_node_predicates
from pybel.struct.filters.typing import NodePredicates

__all__ = [
    'group_nodes_by_annotation',
    'average_node_annotation',
    'group_nodes_by_annotation_filtered'
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
        node_predicates: NodePredicates = None,
        annotation: str = 'Subgraph',
) -> Mapping[str, Set[BaseEntity]]:
    """Group the nodes occurring in edges by the given annotation, with a node filter applied.

    :param graph: A BEL graph
    :param node_predicates: A predicate or list of predicates (graph, node) -> bool
    :param annotation: The annotation to use for grouping
    :return: A dictionary of {annotation value: set of nodes}
    """
    node_filter = concatenate_node_predicates(node_predicates)

    return {
        key: {
            node
            for node in nodes
            if node_filter(graph, node)
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
    parent_to_variants = defaultdict(set)
    names = set(names)

    for u, v, d in graph.edges(data=True):
        if d[RELATION] in {HAS_MEMBER, HAS_COMPONENT} and v.get(NAMESPACE) == namespace and v.get(NAME) in names:
            parent_to_variants[v].add(u)

        elif d[RELATION] == HAS_VARIANT and u.get(NAMESPACE) == namespace and u.get(NAME) in names:
            parent_to_variants[u].add(v)

        elif d[RELATION] == ORTHOLOGOUS and u.get(NAMESPACE) == namespace and u.get(NAME) in names:
            parent_to_variants[u].add(v)

    return dict(parent_to_variants)
