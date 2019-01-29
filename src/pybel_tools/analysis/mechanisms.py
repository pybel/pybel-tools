# -*- coding: utf-8 -*-

"""This module contains functions to compare generated subgraphs to canonical subgraphs."""

import itertools as itt
from collections import defaultdict
from typing import Any, Dict, Mapping

from pybel import BELGraph
from pybel.struct import get_subgraphs_by_annotation
from ..generation import generate_bioprocess_mechanisms
from ..utils import tanimoto_set_similarity

__all__ = [
    'compare'
]


def compare(graph: BELGraph, annotation: str = 'Subgraph') -> Mapping[str, Mapping[str, float]]:
    """Compare generated mechanisms to actual ones.

    1. Generates candidate mechanisms for each biological process
    2. Gets sub-graphs for all NeuroMMSig signatures
    3. Make tanimoto similarity comparison for all sets

    :return: A dictionary table comparing the canonical subgraphs to generated ones
    """
    canonical_mechanisms = get_subgraphs_by_annotation(graph, annotation)
    canonical_nodes = _transform_graph_dict_to_node_dict(canonical_mechanisms)

    candidate_mechanisms = generate_bioprocess_mechanisms(graph)
    candidate_nodes = _transform_graph_dict_to_node_dict(candidate_mechanisms)

    results: Dict[str, Dict[str, float]] = defaultdict(dict)

    it = itt.product(canonical_nodes.items(), candidate_nodes.items())
    for (canonical_name, canonical_graph), (candidate_bp, candidate_graph) in it:
        tanimoto = tanimoto_set_similarity(candidate_nodes, canonical_nodes)
        results[canonical_name][candidate_bp] = tanimoto

    return dict(results)


def _transform_graph_dict_to_node_dict(d: Mapping[Any, BELGraph]):
    return {key: set(graph) for key, graph in d.items()}
