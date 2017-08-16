# -*- coding: utf-8 -*-

"""This module contains functions that help add description and node labels"""

from . import go_annotator
from . import hgnc
from . import node_annotator
from .go_annotator import *
from .hgnc import *
from .node_annotator import *

__all__ = (
    node_annotator.__all__ +
    hgnc.__all__ +
    go_annotator.__all__
)
