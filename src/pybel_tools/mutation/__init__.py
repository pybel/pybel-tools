# -*- coding: utf-8 -*-

"""Mutation functions to supplement :mod:`pybel.struct.mutation`."""

from . import collapse, deletion, expansion, highlight, inference, metadata, random
from .bound import *
from .collapse import *
from .deletion import *
from .expansion import *
from .highlight import *
from .inference import *
from .metadata import *
from .random import *

__all__ = (
    collapse.__all__ +
    deletion.__all__ +
    expansion.__all__ +
    highlight.__all__ +
    inference.__all__ +
    metadata.__all__ +
    random.__all__ +
    bound.__all__
)
