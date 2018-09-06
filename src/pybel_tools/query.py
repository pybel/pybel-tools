# -*- coding: utf-8 -*-

"""A wrapper around the PyBEL query builder."""

import warnings

from pybel.struct.query import Query, QueryMissingNetworksError, SEED_DATA, SEED_METHOD

__all__ = [
    'QueryMissingNetworksError',
    'Query',
    'SEED_DATA',
    'SEED_METHOD',
]

warnings.warn('Query functionality has moved to pybel.struct.query', DeprecationWarning)
