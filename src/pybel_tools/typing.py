# -*- coding: utf-8 -*-

"""Type hints for PyBEL-Tools."""

from typing import Set, Tuple

from pybel import BaseEntity

__all__ = [
    'EdgeSet',
    'NodePair',
    'NodeTriple',
    'SetOfNodePairs',
    'SetOfNodeTriples',
]

NodePair = Tuple[BaseEntity, BaseEntity]
EdgeSet = SetOfNodePairs = Set[NodePair]

NodeTriple = Tuple[BaseEntity, BaseEntity, BaseEntity]
SetOfNodeTriples = Set[NodeTriple]
