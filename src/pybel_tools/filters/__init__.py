# -*- coding: utf-8 -*-

"""

This module contains functions for filtering node and edge iterables. It relies heavily on the concepts of
`functional programming <https://en.wikipedia.org/wiki/Functional_programming>`_ and the concept of
`predicates <https://en.wikipedia.org/wiki/Predicate_(mathematical_logic)>`_.

"""

from . import edge_filters
from . import node_deletion
from . import node_filters
from . import node_selection
from .edge_filters import *
from .node_deletion import *
from .node_filters import *
from .node_selection import *

__all__ = (
    edge_filters.__all__ +
    node_deletion.__all__ +
    node_filters.__all__ +
    node_selection.__all__
)
