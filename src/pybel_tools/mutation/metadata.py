# -*- coding: utf-8 -*-

"""Metadata mutation utilities."""

import logging
from typing import Set

from pybel import BELGraph, Manager
from pybel.constants import CITATION, IDENTIFIER
from pybel.manager.citation_utils import enrich_pubmed_citations, get_citations_by_pmids
from pybel.struct.filters import filter_edges, has_pubmed
from pybel.struct.pipeline import in_place_transformation
from pybel.struct.summary import get_pubmed_identifiers

__all__ = [
    'enrich_pubmed_citations',
]

logger = logging.getLogger(__name__)


@in_place_transformation
def enrich_pubmed_citations(graph: BELGraph, manager: Manager) -> Set[str]:
    """Overwrite all PubMed citations with values from NCBI's eUtils lookup service.

    :return: A set of PMIDs for which the eUtils service crashed
    """
    pmids = get_pubmed_identifiers(graph)
    pmid_data, errors = get_citations_by_pmids(manager=manager, pmids=pmids)

    for u, v, k in filter_edges(graph, has_pubmed):
        pmid = graph[u][v][k][CITATION].identifier

        if pmid not in pmid_data:
            logger.warning('Missing data for PubMed identifier: %s', pmid)
            errors.add(pmid)
            continue

        graph[u][v][k][CITATION].update(pmid_data[pmid])

    return errors
