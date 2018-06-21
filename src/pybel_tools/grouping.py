# -*- coding: utf-8 -*-

from collections import defaultdict

import logging

from pybel.constants import ANNOTATIONS
from pybel.struct.utils import update_metadata
from .utils import safe_add_edge

log = logging.getLogger(__name__)

__all__ = [
    'get_subgraphs_by_annotation_filtered',
]


def get_subgraphs_by_annotation_filtered(graph, annotation, values):
    """Stratifies the given graph into subgraphs based on the values for edges' annotations, but filter by a set
    of given values

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to group by
    :param iter[str] values: The values to keep
    :rtype: dict[str,pybel.BELGraph]
    """
    result = defaultdict(graph.fresh_copy)
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
