# -*- coding: utf-8 -*-

"""This module contains functions to summarize the provenance (citations, evidences, and authors) in a BEL graph"""

import itertools as itt
import logging
from collections import Counter, defaultdict
from datetime import datetime

from pybel.constants import *
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


def _generate_citation_dict(graph):
    """Prepares a citation data dictionary from a graph. This is a helper function

    :param pybel.BELGraph graph: A BEL graph
    :return: A dictionary of {citation type: {(reference, name) -> {set of (source node, target node)}}}
    :rtype: dict[str,dict[tuple[tuple,tuple],str]]
    """
    results = defaultdict(lambda: defaultdict(set))

    for u, v, data in graph.edges_iter(data=True):
        if CITATION not in data:
            continue
        results[data[CITATION][CITATION_TYPE]][u, v].add(data[CITATION][CITATION_REFERENCE].strip())

    return dict(results)


def get_pmid_by_keyword(keyword, graph=None, pubmed_identifiers=None):
    """Gets the set of PubMed identifiers beginning with the given keyword string
    
    :param pybel.BELGraph graph: A BEL graph
    :param str keyword: The beginning of a PubMed identifier
    :param set[str] pubmed_identifiers: A set of pre-cached PubMed identifiers
    :return: A set of PubMed identifiers starting with the given string
    :rtype: set[str]
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


def count_pmids(graph):
    """Counts the frequency of PubMed documents in a graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A Counter from {(pmid, name): frequency}
    :rtype: collections.Counter
    """
    return Counter(iterate_pubmed_identifiers(graph))


def get_citation_pair(data):
    return data[CITATION][CITATION_TYPE], data[CITATION][CITATION_REFERENCE].strip()


def count_unique_citations(graph):
    """Returns the number of unique citations

    :param pybel.BELGraph graph: A BEL graph
    :return: The number of unique citations in the graph.
    :rtype: int
    """
    return len({
        get_citation_pair(data)
        for data in graph_edge_data_iter(graph)
        if CITATION in data
    })


def count_citations(graph, **annotations):
    """Counts the citations in a graph based on a given filter

    :param pybel.BELGraph graph: A BEL graph
    :param dict annotations: The annotation filters to use
    :return: A counter from {(citation type, citation reference): frequency}
    :rtype: collections.Counter
    """
    citations = defaultdict(set)

    annotation_dict_filter = build_edge_data_filter(annotations)

    for u, v, _, d in filter_edges(graph, annotation_dict_filter):
        if CITATION not in d:
            continue

        citations[u, v].add(get_citation_pair(d))

    counter = Counter(itt.chain.from_iterable(citations.values()))
    return counter


def count_citations_by_annotation(graph, annotation):
    """Groups the citation counters by subgraphs induced by the annotation

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to use to group the graph
    :return: A dictionary of Counters {subgraph name: Counter from {citation: frequency}}
    """
    citations = defaultdict(lambda: defaultdict(set))
    for u, v, data in graph.edges_iter(data=True):
        if not edge_has_annotation(data, annotation) or CITATION not in data:
            continue

        k = data[ANNOTATIONS][annotation]

        citations[k][u, v].add((data[CITATION][CITATION_TYPE], data[CITATION][CITATION_REFERENCE].strip()))

    return {k: Counter(itt.chain.from_iterable(v.values())) for k, v in citations.items()}


def check_authors_in_data(data):
    return CITATION not in data or CITATION_AUTHORS not in data[CITATION]


def raise_for_unparsed_authors(data):
    authors = data[CITATION][CITATION_AUTHORS]
    if isinstance(authors, str):
        raise ValueError('Graph should be converted with pbt.mutation.parse_authors first: {}'.format(authors))


def count_authors(graph):
    """Counts the contributions of each author to the given graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A Counter from {author name: frequency}
    :rtype: collections.Counter
    """
    authors = []
    for data in graph_edge_data_iter(graph):
        if check_authors_in_data(data):
            continue
        raise_for_unparsed_authors(data)
        for author in data[CITATION][CITATION_AUTHORS]:
            authors.append(author)

    return Counter(authors)


def count_author_publications(graph):
    """Counts the number of publications of each author to the given graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A Counter from {author name: frequency}
    :rtype: collections.Counter
    """
    authors = defaultdict(list)
    for data in graph_edge_data_iter(graph):
        if check_authors_in_data(data):
            continue
        raise_for_unparsed_authors(data)
        for author in data[CITATION][CITATION_AUTHORS]:
            authors[author].append(data[CITATION][CITATION_REFERENCE].strip())

    return Counter(count_dict_values(count_defaultdict(authors)))


# TODO switch to use node filters
def get_authors(graph):
    """Gets the set of all authors in the given graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A set of author names
    :rtype: set[str]
    """
    result = set()

    for data in graph_edge_data_iter(graph):
        if check_authors_in_data(data):
            continue

        authors = data[CITATION][CITATION_AUTHORS]

        result.update(
            authors.strip().split('|')
            if isinstance(authors, str)
            else authors
        )

    return result


def count_unique_authors(graph):
    """Counts all authors in the given graph

    :param pybel.BELGraph graph: A BEL graph
    :return: The number of unique authors whose publications contributed to the graph
    :rtype: int
    """
    return len(get_authors(graph))


def get_authors_by_keyword(keyword, graph=None, authors=None):
    """Gets authors for whom the search term is a substring
    
    :param pybel.BELGraph graph: A BEL graph
    :param str keyword: The keyword to search the author strings for
    :param set[str] authors: An optional set of pre-cached authors calculated from the graph
    :return: A set of authors with the keyword as a substring
    :rtype: set[str]
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


def count_authors_by_annotation(graph, annotation='Subgraph'):
    """Groups the author counters by subgraphs induced by the annotation

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to use to group the graph
    :return: A dictionary of Counters {subgraph name: Counter from {author: frequency}}
    :rtype: dict
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


def get_evidences_by_pmid(graph, pmids):
    """Gets a dictionary from the given PubMed identifiers to the sets of all evidence strings associated with each
    in the graph

    :param pybel.BELGraph graph: A BEL graph
    :param str or iter[str] pmids: An iterable of PubMed identifiers, as strings. Is consumed and converted to a set.
    :return: A dictionary of {pmid: set of all evidence strings}
    :rtype: dict
    """
    result = defaultdict(set)

    for _, _, _, data in filter_edges(graph, build_pmid_inclusion_filter(pmids)):
        result[data[CITATION][CITATION_REFERENCE]].add(data[EVIDENCE])

    return dict(result)


# TODO date parsing should be handled during either pybel parse-time or during graph loading.
def count_citation_years(graph):
    """Counts the number of citations in each year

    :param pybel.BELGraph graph: A BEL graph
    :return: A Counter of {int year: int frequency}
    :rtype: collections.Counter
    """
    result = defaultdict(set)

    for data in graph_edge_data_iter(graph):
        if CITATION not in data or CITATION_DATE not in data[CITATION]:
            continue

        try:
            dt = _ensure_datetime(data[CITATION][CITATION_DATE])
            result[dt.year].add((data[CITATION][CITATION_TYPE], data[CITATION][CITATION_REFERENCE]))
        except:
            continue

    return count_dict_values(result)


def _ensure_datetime(s):
    if isinstance(s, datetime):
        return s

    elif isinstance(s, str):
        return datetime.strptime(s, '%Y-%m-%d')

    raise TypeError


def create_timeline(year_counter):
    """Completes the Counter timeline

    :param Counter year_counter: counter dict for each year
    :return: complete timeline
    :rtype: list[tuple[int,int]]
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


def get_citation_years(graph):
    """Creates a citation timeline counter

    :param pybel.BELGraph graph: A BEL graph
    :rtype: list[tuple[int,int]]
    """
    return create_timeline(count_citation_years(graph))
