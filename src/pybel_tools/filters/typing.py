# -*- coding: utf-8 -*-

"""Types for filters."""

from typing import Callable, Iterable, Union

from pybel import BELGraph
from pybel.dsl import BaseEntity

__all__ = [
    'NodePredicate',
    'NodePredicates',
    'EdgePredicate',
    'EdgePredicates',
    'Strings',
]

NodePredicate = Callable[[BELGraph, BaseEntity], bool]
NodePredicates = Union[None, NodePredicate, Iterable[NodePredicate]]

EdgePredicate = Callable[[BELGraph, BaseEntity, BaseEntity, str], bool]
EdgePredicates = Union[None, EdgePredicate, Iterable[EdgePredicate]]

Strings = Union[str, Iterable[str]]
