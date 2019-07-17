# -*- coding: utf-8 -*-

"""This module builds mutation functions that are bound to a manager."""

from typing import Callable

from pybel import BELGraph, Manager
from pybel.struct.mutation.expansion.neighborhood import expand_node_neighborhood
from pybel.struct.pipeline import in_place_transformation, uni_in_place_transformation

__all__ = [
    'build_expand_node_neighborhood_by_hash',
    'build_delete_node_by_hash',
]


def build_expand_node_neighborhood_by_hash(manager: Manager) -> Callable[[BELGraph, BELGraph, str], None]:  # noqa: D202
    """Make an expand function that's bound to the manager."""

    @uni_in_place_transformation
    def expand_node_neighborhood_by_hash(universe: BELGraph, graph: BELGraph, node_hash: str) -> None:
        """Expand around the neighborhoods of a node by identifier."""
        node = manager.get_dsl_by_hash(node_hash)
        return expand_node_neighborhood(universe, graph, node)

    return expand_node_neighborhood_by_hash


def build_delete_node_by_hash(manager: Manager) -> Callable[[BELGraph, str], None]:  # noqa: D202
    """Make a delete function that's bound to the manager."""

    @in_place_transformation
    def delete_node_by_hash(graph: BELGraph, node_hash: str) -> None:
        """Remove a node by identifier."""
        node = manager.get_dsl_by_hash(node_hash)
        graph.remove_node(node)

    return delete_node_by_hash
