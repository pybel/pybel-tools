# -*- coding: utf-8 -*-

"""This module contains functions to summarize the provenance (citations, evidences, and authors) in a BEL graph."""

import itertools as itt
import logging
import typing
from collections import Counter, defaultdict
from datetime import datetime
from typing import Iterable, List, Mapping, Optional, Set, Tuple, Union

from pybel import BELGraph
from pybel.constants import (
    ANNOTATIONS, CITATION, CITATION_AUTHORS, CITATION_DATE, CITATION_REFERENCE, CITATION_TYPE, EVIDENCE,
)
from pybel.dsl import BaseEntity
from pybel.struct.filters import filter_edges
from pybel.struct.filters.edge_predicates import edge_has_annotation
from pybel.struct.summary import iterate_pubmed_identifiers
from pybel.typing import Strings
from ..filters import build_edge_data_filter, build_pmid_inclusion_filter
from ..utils import count_defaultdict, count_dict_values, group_as_lists, group_as_sets

__all__ = [
    'count_pmids',
    'get_pmid_by_keyword',
    'count_citations',
    'count_citations_by_annotation',
    'count_authors',
    'count_author_publications',
    'get_authors',
    'get_authors_by_keyword',
    'count_authors_by_annotation',
    'get_evidences_by_pmid',
    'count_citation_years',
    'create_timeline',
    'get_citation_years',
    'count_confidences',
]

log = logging.getLogger(__name__)


def _generate_citation_dict(graph: BELGraph) -> Mapping[str, Mapping[Tuple[BaseEntity, BaseEntity], str]]:
    """Prepare a citation data dictionary from a graph.

    :return: A dictionary of dictionaries {citation type: {(source, target): citation reference}
    """
    results = defaultdict(lambda: defaultdict(set))

    for u, v, data in graph.edges(data=True):
        if CITATION not in data:
            continue
        results[data[CITATION][CITATION_TYPE]][u, v].add(data[CITATION][CITATION_REFERENCE].strip())

    return dict(results)


def get_pmid_by_keyword(
        keyword: str,
        graph: Optional[BELGraph] = None,
        pubmed_identifiers: Optional[Set[str]] = None,
) -> Set[str]:
    """Get the set of PubMed identifiers beginning with the given keyword string.

    :param keyword: The beginning of a PubMed identifier
    :param graph: A BEL graph
    :param pubmed_identifiers: A set of pre-cached PubMed identifiers
    :return: A set of PubMed identifiers starting with the given string
    """
    if pubmed_identifiers is not None:
        return {
            pubmed_identifier
            for pubmed_identifier in pubmed_identifiers
            if pubmed_identifier.startswith(keyword)
        }

    if graph is None:
        raise ValueError('Graph not supplied')

    return {
        pubmed_identifier
        for pubmed_identifier in iterate_pubmed_identifiers(graph)
        if pubmed_identifier.startswith(keyword)
    }


def count_pmids(graph: BELGraph) -> Counter:
    """Count the frequency of PubMed documents in a graph.

    :return: A Counter from {(pmid, name): frequency}
    """
    return Counter(iterate_pubmed_identifiers(graph))


def count_citations(graph: BELGraph, **annotations) -> Counter:
    """Count the citations in a graph based on a given filter.

    :param graph: A BEL graph
    :param dict annotations: The annotation filters to use
    :return: A counter from {(citation type, citation reference): frequency}
    """
    annotation_dict_filter = build_edge_data_filter(annotations)

    citations = defaultdict(set)
    for u, v, _, d in filter_edges(graph, annotation_dict_filter):
        if CITATION in d:
            citations[u, v].add((d[CITATION][CITATION_TYPE], d[CITATION][CITATION_REFERENCE].strip()))

    return Counter(itt.chain.from_iterable(citations.values()))


def count_citations_by_annotation(graph: BELGraph, annotation: str) -> Mapping[str, typing.Counter[str]]:
    """Group the citation counters by subgraphs induced by the annotation.

    :param graph: A BEL graph
    :param annotation: The annotation to use to group the graph
    :return: A dictionary of Counters {subgraph name: Counter from {citation: frequency}}
    """
    citations = defaultdict(lambda: defaultdict(set))
    for u, v, data in graph.edges(data=True):
        if not edge_has_annotation(data, annotation) or CITATION not in data:
            continue

        k = data[ANNOTATIONS][annotation]

        citations[k][u, v].add((data[CITATION][CITATION_TYPE], data[CITATION][CITATION_REFERENCE].strip()))

    return {
        k: Counter(itt.chain.from_iterable(v.values()))
        for k, v in citations.items()
    }


def count_authors(graph: BELGraph) -> typing.Counter[str]:
    """Count the number of edges in which each author appears."""
    return Counter(graph._iterate_authors())


def count_author_publications(graph: BELGraph) -> typing.Counter[str]:
    """Count the number of publications of each author to the given graph."""
    authors = group_as_lists(_iter_author_publiations(graph))
    return Counter(count_dict_values(count_defaultdict(authors)))


def _iter_author_publiations(graph: BELGraph) -> Iterable[Tuple[str, Tuple[str, str]]]:
    for _, _, data in graph.edges(data=True):
        if CITATION not in data or CITATION_AUTHORS not in data[CITATION]:
            continue
        for author in data[CITATION][CITATION_AUTHORS]:
            yield (
                author,
                (data[CITATION][CITATION_TYPE], data[CITATION][CITATION_REFERENCE])
            )


def get_authors(graph: BELGraph) -> Set[str]:
    """Get the set of all authors in the given graph."""
    return set(graph._iterate_authors())


def get_authors_by_keyword(keyword: str, graph=None, authors=None) -> Set[str]:
    """Get authors for whom the search term is a substring.

    :param pybel.BELGraph graph: A BEL graph
    :param keyword: The keyword to search the author strings for
    :param set[str] authors: An optional set of pre-cached authors calculated from the graph
    :return: A set of authors with the keyword as a substring
    """
    keyword_lower = keyword.lower()

    if authors is not None:
        return {
            author
            for author in authors
            if keyword_lower in author.lower()
        }

    if graph is None:
        raise ValueError('Graph not supplied')

    return {
        author
        for author in get_authors(graph)
        if keyword_lower in author.lower()
    }


def count_authors_by_annotation(graph: BELGraph, annotation: str = 'Subgraph') -> Mapping[str, typing.Counter[str]]:
    """Group the author counters by sub-graphs induced by the annotation.

    :param graph: A BEL graph
    :param annotation: The annotation to use to group the graph
    :return: A dictionary of Counters {subgraph name: Counter from {author: frequency}}
    """
    authors = group_as_lists(_iter_authors_by_annotation(graph, annotation=annotation))
    return count_defaultdict(authors)


def _iter_authors_by_annotation(graph: BELGraph, annotation: str = 'Subgraph') -> Iterable[Tuple[str, str]]:
    for _, _, data in graph.edges(data=True):
        if not edge_has_annotation(data, annotation) or CITATION not in data or CITATION_AUTHORS not in data[CITATION]:
            continue
        for author in data[CITATION][CITATION_AUTHORS]:
            yield data[ANNOTATIONS][annotation], author


def get_evidences_by_pmid(graph: BELGraph, pmids: Strings) -> Mapping[str, Set[str]]:
    """Map PubMed identifiers to their evidence strings appearing in the graph.

    :param graph: A BEL graph
    :param pmids: An iterable of PubMed identifiers, as strings. Is consumed and converted to a set.
    :return: A dictionary of {pmid: set of all evidence strings}
    :rtype: dict
    """
    return group_as_sets(
        (data[CITATION][CITATION_REFERENCE], data[EVIDENCE])
        for _, _, _, data in filter_edges(graph, build_pmid_inclusion_filter(pmids))
    )


# TODO date parsing should be handled during either pybel parse-time or during graph loading.
def count_citation_years(graph: BELGraph) -> typing.Counter[int]:
    """Count the number of citations from each year."""
    result = defaultdict(set)

    for _, _, data in graph.edges(data=True):
        if CITATION not in data or CITATION_DATE not in data[CITATION]:
            continue

        try:
            dt = _ensure_datetime(data[CITATION][CITATION_DATE])
        except ValueError:
            continue
        else:
            result[dt.year].add((data[CITATION][CITATION_TYPE], data[CITATION][CITATION_REFERENCE]))

    return count_dict_values(result)


def _ensure_datetime(s: Union[datetime, str]) -> datetime:
    if isinstance(s, datetime):
        return s

    elif isinstance(s, str):
        return datetime.strptime(s, '%Y-%m-%d')

    raise TypeError


def get_citation_years(graph: BELGraph) -> List[Tuple[int, int]]:
    """Create a citation timeline counter from the graph."""
    return create_timeline(count_citation_years(graph))


def create_timeline(year_counter: typing.Counter[int]) -> List[Tuple[int, int]]:
    """Complete the Counter timeline.

    :param Counter year_counter: counter dict for each year
    :return: complete timeline
    """
    if not year_counter:
        return []

    from_year = min(year_counter) - 1
    until_year = datetime.now().year + 1

    return [
        (year, year_counter.get(year, 0))
        for year in range(from_year, until_year)
    ]


def count_confidences(graph: BELGraph) -> typing.Counter[str]:
    """Count the confidences in the graph."""
    return Counter(
        (
            'None'
            if ANNOTATIONS not in data or 'Confidence' not in data[ANNOTATIONS] else
            list(data[ANNOTATIONS]['Confidence'])[0]
        )
        for _, _, data in graph.edges(data=True)
        if CITATION in data  # don't bother with unqualified statements
    )
