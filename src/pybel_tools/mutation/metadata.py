# -*- coding: utf-8 -*-

import logging

from typing import Set

from pybel.constants import CITATION, CITATION_REFERENCE
from pybel.manager.citation_utils import get_citations_by_pmids
from pybel.struct.filters import filter_edges
from pybel.struct.filters.edge_predicates import has_pubmed
from pybel.struct.pipeline import in_place_transformation, uni_in_place_transformation
from pybel.struct.summary import get_pubmed_identifiers
from pybel.struct.summary.node_summary import get_namespaces
from ..summary.edge_summary import get_annotations

__all__ = [
    'enrich_pubmed_citations',
]

log = logging.getLogger(__name__)


@in_place_transformation
def enrich_pubmed_citations(graph, stringify_authors:bool=False, manager=None) -> Set[str]:
    """Overwrites all PubMed citations with values from NCBI's eUtils lookup service.

    Sets authors as list, so probably a good idea to run :func:`pybel_tools.mutation.serialize_authors` before
    exporting.

    :param pybel.BELGraph graph: A BEL graph
    :param bool stringify_authors: Converts all author lists to author strings using
                                  :func:`pybel_tools.mutation.serialize_authors`. Defaults to ``False``.
    :param manager: An RFC-1738 database connection string, a pre-built :class:`pybel.manager.Manager`,
                    or ``None`` for default connection
    :type manager: None or str or Manager
    :return: A set of PMIDs for which the eUtils service crashed
    """
    if 'PYBEL_ENRICHED_CITATIONS' in graph.graph:
        log.warning('citations have already been enriched in %s', graph)
        return set()

    pmids = get_pubmed_identifiers(graph)
    pmid_data, errors = get_citations_by_pmids(manager=manager, pmids=pmids)

    for u, v, k in filter_edges(graph, has_pubmed):
        pmid = graph[u][v][k][CITATION][CITATION_REFERENCE].strip()

        if pmid not in pmid_data:
            log.warning('Missing data for PubMed identifier: %s', pmid)
            errors.add(pmid)
            continue

        graph[u][v][k][CITATION].update(pmid_data[pmid])

    graph.graph['PYBEL_ENRICHED_CITATIONS'] = True

    return errors


@uni_in_place_transformation
def update_context(universe, graph):
    """Updates the context of a subgraph from the universe of all knowledge.

    :param pybel.BELGraph universe: The universe of knowledge
    :param pybel.BELGraph graph: A BEL graph
    """
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
