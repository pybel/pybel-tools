# -*- coding: utf-8 -*-

import logging
from collections import defaultdict

from pybel import BELGraph
from pybel.constants import ANNOTATIONS
from . import pipeline
from .selection.induce_subgraph import update_metadata
from .utils import safe_add_edge

log = logging.getLogger(__name__)

__all__ = [
    'get_subgraphs_by_annotation',
    'get_subgraphs_by_annotation_filtered',
]


@pipeline.splitter
def get_subgraphs_by_annotation(graph, annotation, keep_undefined=True, sentinel='Undefined'):
    """Stratifies the given graph into subgraphs based on the values for edges' annotations

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to group by
    :param bool keep_undefined: If true, uses the sentinel value to store a subgraph of edges not matching the given
     annotation.
    :param str sentinel: The value to stick unannotated edges into
    :rtype: dict[str,pybel.BELGraph]
    """
    result = defaultdict(BELGraph)

    for source, target, key, data in graph.edges_iter(keys=True, data=True):
        annotation_dict = data.get(ANNOTATIONS)

        if annotation_dict is None or annotation not in annotation_dict:
            if keep_undefined:
                safe_add_edge(result[sentinel], source, target, key, data)
        else:
            for value in annotation_dict[annotation]:
                safe_add_edge(result[value], source, target, key, data)

    for value in result.values():
        update_metadata(value, graph)

    return dict(result)


@pipeline.splitter
def get_subgraphs_by_annotation_filtered(graph, annotation, values):
    """Stratifies the given graph into subgraphs based on the values for edges' annotations, but filter by a set
    of given values

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to group by
    :param iter[str] values: The values to keep
    :rtype: dict[str,pybel.BELGraph]
    """
    result = defaultdict(BELGraph)
    values = set(values)

    for source, target, key, data in graph.edges_iter(keys=True, data=True):
        annotation_dict = data.get(ANNOTATIONS)

        if annotation_dict is None or annotation not in annotation_dict:
            continue

        for value in annotation_dict[annotation]:
            if value in values:
                safe_add_edge(result[value], source, target, key, data)

    for value in result.values():
        update_metadata(value, graph)

    return dict(result)
