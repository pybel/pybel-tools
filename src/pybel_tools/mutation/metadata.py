# -*- coding: utf-8 -*-

"""Metadata mutation utilities."""

import logging
from typing import Set

from pybel import BELGraph, Manager
from pybel.constants import CITATION, CITATION_REFERENCE
from pybel.manager.citation_utils import get_citations_by_pmids
from pybel.struct.filters import filter_edges, has_pubmed
from pybel.struct.pipeline import in_place_transformation, uni_in_place_transformation
from pybel.struct.summary import get_annotations, get_namespaces, get_pubmed_identifiers

__all__ = [
    'enrich_pubmed_citations',
]

log = logging.getLogger(__name__)


@in_place_transformation
def enrich_pubmed_citations(graph: BELGraph, manager: Manager) -> Set[str]:
    """Overwrite all PubMed citations with values from NCBI's eUtils lookup service.

    :return: A set of PMIDs for which the eUtils service crashed
    """
    pmids = get_pubmed_identifiers(graph)
    pmid_data, errors = get_citations_by_pmids(manager=manager, pmids=pmids)

    for u, v, k in filter_edges(graph, has_pubmed):
        pmid = graph[u][v][k][CITATION][CITATION_REFERENCE].strip()

        if pmid not in pmid_data:
            log.warning('Missing data for PubMed identifier: %s', pmid)
            errors.add(pmid)
            continue

        graph[u][v][k][CITATION].update(pmid_data[pmid])

    return errors


@uni_in_place_transformation
def update_context(universe: BELGraph, graph: BELGraph):
    """Update the context of a subgraph from the universe of all knowledge."""
    for namespace in get_namespaces(graph):
        if namespace in universe.namespace_url:
            graph.namespace_url[namespace] = universe.namespace_url[namespace]
        elif namespace in universe.namespace_pattern:
            graph.namespace_pattern[namespace] = universe.namespace_pattern[namespace]
        else:
            log.warning('namespace: %s missing from universe', namespace)

    for annotation in get_annotations(graph):
        if annotation in universe.annotation_url:
            graph.annotation_url[annotation] = universe.annotation_url[annotation]
        elif annotation in universe.annotation_pattern:
            graph.annotation_pattern[annotation] = universe.annotation_pattern[annotation]
        elif annotation in universe.annotation_list:
            graph.annotation_list[annotation] = universe.annotation_list[annotation]
        else:
            log.warning('annotation: %s missing from universe', annotation)
