# -*- coding: utf-8 -*-

from .deletion import prune_central_dogma
from .inference import infer_central_dogma
from .. import pipeline


@pipeline.in_place_mutator
def opening_on_central_dogma(graph):
    """Infers central dogmatic relations with :func:`infer_central_dogma` then successively prunes gene leaves then
    RNA leaves with :func:`prune_central_dogma` to connect disparate elements in a knowledge assembly

    :param pybel.BELGraph graph: A BEL graph

    Equivalent to:

    >>> infer_central_dogma(graph)
    >>> prune_central_dogma(graph)
    """
    infer_central_dogma(graph)
    prune_central_dogma(graph)
