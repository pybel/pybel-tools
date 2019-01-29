# -*- coding: utf-8 -*-

"""This module contains functions to summarize the provenance (citations, evidences, and authors) in a BEL graph."""

import logging
from collections import Counter, defaultdict
from datetime import datetime
from typing import Iterable, List, Mapping, Optional, Set, Tuple, Union

import itertools as itt

from pybel import BELGraph
from pybel.constants import *
from pybel.dsl import BaseEntity
from pybel.struct.filters import filter_edges
from pybel.struct.filters.edge_predicates import edge_has_annotation
from pybel.struct.summary import iterate_pubmed_identifiers
from ..filters import build_edge_data_filter, build_pmid_inclusion_filter
from ..utils import count_defaultdict, count_dict_values, graph_edge_data_iter

__all__ = [
    'count_pmids',
    'get_pmid_by_keyword',
    'count_citations',
    'count_citations_by_annotation',
    'count_authors',
    'count_unique_authors',
    'count_author_publications',
    'count_unique_citations',
    'get_authors',
    'get_authors_by_keyword',
    'count_authors_by_annotation',
    'get_evidences_by_pmid',
    'count_citation_years',
    'create_timeline',
    'get_citation_years'
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


def get_pmid_by_keyword(keyword: str,
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


def get_citation_pair(data):
    return data[CITATION][CITATION_TYPE], data[CITATION][CITATION_REFERENCE].strip()


def count_unique_citations(graph: BELGraph) -> int:
    """Return the number of unique citations in the graph."""
    return len({
        get_citation_pair(data)
        for data in graph_edge_data_iter(graph)
        if CITATION in data
    })


def count_citations(graph: BELGraph, **annotations) -> Counter:
    """Counts the citations in a graph based on a given filter

    :param graph: A BEL graph
    :param dict annotations: The annotation filters to use
    :return: A counter from {(citation type, citation reference): frequency}
    """
    citations = defaultdict(set)

    annotation_dict_filter = build_edge_data_filter(annotations)

    for u, v, _, d in filter_edges(graph, annotation_dict_filter):
        if CITATION not in d:
            continue

        citations[u, v].add(get_citation_pair(d))

    return Counter(itt.chain.from_iterable(citations.values()))


def count_citations_by_annotation(graph: BELGraph, annotation: str) -> Mapping[str, Counter]:
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

    return {k: Counter(itt.chain.from_iterable(v.values())) for k, v in citations.items()}


def _check_authors_in_data(data) -> bool:
    return CITATION not in data or CITATION_AUTHORS not in data[CITATION]


def _raise_for_unparsed_authors(data: Mapping) -> None:
    authors = data[CITATION][CITATION_AUTHORS]
    if isinstance(authors, str):
        raise ValueError('Graph should be converted with pbt.mutation.parse_authors first: {}'.format(authors))


def count_authors(graph: BELGraph) -> Counter:
    """Count the number of edges in which each author appears."""
    authors = []

    for data in graph_edge_data_iter(graph):
        if _check_authors_in_data(data):
            continue
        _raise_for_unparsed_authors(data)
        for author in data[CITATION][CITATION_AUTHORS]:
            authors.append(author)

    return Counter(authors)


def count_author_publications(graph: BELGraph) -> Counter:
    """Count the number of publications of each author to the given graph."""
    authors = defaultdict(list)

    for data in graph_edge_data_iter(graph):
        if _check_authors_in_data(data):
            continue
        _raise_for_unparsed_authors(data)
        for author in data[CITATION][CITATION_AUTHORS]:
            authors[author].append(data[CITATION][CITATION_REFERENCE].strip())

    return Counter(count_dict_values(count_defaultdict(authors)))


# TODO switch to use node filters
def get_authors(graph: BELGraph) -> Set[str]:
    """Get the set of all authors in the given graph."""
    return set(_iterate_authors(graph))


def _iterate_authors(graph: BELGraph) -> Iterable[str]:
    for _, _, data in graph.edges(data=True):
        if _check_authors_in_data(data):
            continue

        authors = data[CITATION][CITATION_AUTHORS]

        if isinstance(authors, str):
            yield authors
        else:
            yield from authors.strip().split('|')


def count_unique_authors(graph: BELGraph) -> int:
    """Count the number of unique authors whose publications contributed to the graph."""
    return len(get_authors(graph))


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


def count_authors_by_annotation(graph: BELGraph, annotation: str = 'Subgraph') -> Mapping[str, Counter]:
    """Groups the author counters by subgraphs induced by the annotation

    :param graph: A BEL graph
    :param annotation: The annotation to use to group the graph
    :return: A dictionary of Counters {subgraph name: Counter from {author: frequency}}
    """
    authors = defaultdict(list)

    for data in graph_edge_data_iter(graph):
        if not edge_has_annotation(data, annotation) or CITATION not in data or CITATION_AUTHORS not in data[CITATION]:
            continue
        if isinstance(data[CITATION][CITATION_AUTHORS], str):
            raise ValueError('Graph should be converted with pybel.mutation.parse_authors first')
        for author in data[CITATION][CITATION_AUTHORS]:
            authors[data[ANNOTATIONS][annotation]].append(author)

    return count_defaultdict(authors)


def get_evidences_by_pmid(graph: BELGraph, pmids: Union[str, Iterable[str]]):
    """Ges a dictionary from the given PubMed identifiers to the sets of all evidence strings associated with each
    in the graph

    :param graph: A BEL graph
    :param pmids: An iterable of PubMed identifiers, as strings. Is consumed and converted to a set.
    :return: A dictionary of {pmid: set of all evidence strings}
    :rtype: dict
    """
    result = defaultdict(set)

    for _, _, _, data in filter_edges(graph, build_pmid_inclusion_filter(pmids)):
        result[data[CITATION][CITATION_REFERENCE]].add(data[EVIDENCE])

    return dict(result)


# TODO date parsing should be handled during either pybel parse-time or during graph loading.
def count_citation_years(graph: BELGraph) -> Counter:
    """Count the number of citations from each year."""
    result = defaultdict(set)

    for data in graph_edge_data_iter(graph):
        if CITATION not in data or CITATION_DATE not in data[CITATION]:
            continue

        try:
            dt = _ensure_datetime(data[CITATION][CITATION_DATE])
            result[dt.year].add((data[CITATION][CITATION_TYPE], data[CITATION][CITATION_REFERENCE]))
        except Exception:
            continue

    return count_dict_values(result)


def _ensure_datetime(s):
    if isinstance(s, datetime):
        return s

    elif isinstance(s, str):
        return datetime.strptime(s, '%Y-%m-%d')

    raise TypeError


def get_citation_years(graph: BELGraph) -> List[Tuple[int, int]]:
    """Create a citation timeline counter from the graph."""
    return create_timeline(count_citation_years(graph))


def create_timeline(year_counter: Counter) -> List[Tuple[int, int]]:
    """Complete the Counter timeline.

    :param Counter year_counter: counter dict for each year
    :return: complete timeline
    """
    if not year_counter:
        return []

    until_year = datetime.now().year
    from_year = min(year_counter)

    timeline = [
        (year, year_counter.get(year, 0))
        for year in range(from_year, until_year)
    ]

    return timeline
