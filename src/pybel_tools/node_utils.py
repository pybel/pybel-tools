# -*- coding: utf-8 -*-

"""Utilities for handling nodes."""

import warnings

from pybel.struct.node_utils import (
    flatten_list_abundance, list_abundance_cartesian_expansion,
    reaction_cartesian_expansion,
)

__all__ = [
    'flatten_list_abundance',
    'list_abundance_cartesian_expansion',
    'reaction_cartesian_expansion',
]

warnings.warn('Use pybel.struct.node_utils', DeprecationWarning)
