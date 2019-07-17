# -*- coding: utf-8 -*-

"""Edge filters to supplement :mod:`pybel.struct.filters.edge_filters`."""

from typing import Iterable, Mapping, Set

from pybel import BELGraph
from pybel.constants import CITATION, CITATION_AUTHORS, CITATION_REFERENCE, NAMESPACE
from pybel.dsl import BaseEntity
from pybel.struct.filters import (
    build_annotation_dict_all_filter, build_annotation_dict_any_filter, count_passed_edge_filter,
)
from pybel.struct.filters.edge_predicates import edge_predicate, has_authors, has_pathology_causal, has_pubmed
from pybel.struct.filters.typing import EdgePredicate, EdgePredicates
from pybel.typing import EdgeData, Strings
from pybel.utils import subdict_matches

__all__ = [
    'summarize_edge_filter',
    'build_edge_data_filter',
    'build_annotation_dict_all_filter',
    'build_annotation_dict_any_filter',
    'build_pmid_inclusion_filter',
    'build_pmid_exclusion_filter',
    'build_author_inclusion_filter',
    'has_pathology_causal',
    'build_source_namespace_filter',
    'build_target_namespace_filter',
]


def summarize_edge_filter(graph: BELGraph, edge_predicates: EdgePredicates) -> None:
    """Print a summary of the number of edges passing a given set of filters."""
    passed = count_passed_edge_filter(graph, edge_predicates)
    e = (
        ', '.join(edge_filter.__name__ for edge_filter in edge_predicates)
        if isinstance(edge_predicates, Iterable) else
        edge_predicates.__name__
    )
    print(f'{passed}/{graph.number_of_edges()} edges passed {e}')


def build_edge_data_filter(annotations: Mapping, partial_match: bool = True) -> EdgePredicate:  # noqa: D202
    """Build a filter that keeps edges whose data dictionaries are super-dictionaries to the given dictionary.

    :param annotations: The annotation query dict to match
    :param partial_match: Should the query values be used as partial or exact matches? Defaults to :code:`True`.
    """

    @edge_predicate
    def annotation_dict_filter(data: EdgeData) -> bool:
        """Match edges with the given dictionary as a sub-dictionary."""
        return subdict_matches(data, annotations, partial_match=partial_match)

    return annotation_dict_filter


def build_pmid_inclusion_filter(pmids: Strings) -> EdgePredicate:
    """Pass for edges with citations whose references are one of the given PubMed identifiers.

    :param pmids: A PubMed identifier or list of PubMed identifiers to filter for
    """
    if isinstance(pmids, str):
        @edge_predicate
        def pmid_inclusion_filter(data: EdgeData) -> bool:
            """Pass only for edges with PubMed citations matching the contained PubMed identifier.

            :return: If the edge has a PubMed citation with the contained PubMed identifier
            """
            return has_pubmed(data) and data[CITATION][CITATION_REFERENCE] == pmids

    elif isinstance(pmids, Iterable):
        pmids = set(pmids)

        @edge_predicate
        def pmid_inclusion_filter(data: EdgeData) -> bool:
            """Pass for edges with PubMed citations matching one of the contained PubMed identifiers.

            :param data: The edge data dictionary
            :return: If the edge has a PubMed citation with one of the contained PubMed identifiers
            """
            return has_pubmed(data) and data[CITATION][CITATION_REFERENCE] in pmids

    else:
        raise TypeError

    return pmid_inclusion_filter


def build_pmid_exclusion_filter(pmids: Strings) -> EdgePredicate:
    """Fail for edges with citations whose references are one of the given PubMed identifiers.

    :param pmids: A PubMed identifier or list of PubMed identifiers to filter against
    """
    if isinstance(pmids, str):
        @edge_predicate
        def pmid_exclusion_filter(data: EdgeData) -> bool:
            """Fail for edges with PubMed citations matching the contained PubMed identifier.

            :return: If the edge has a PubMed citation with the contained PubMed identifier
            """
            return has_pubmed(data) and data[CITATION][CITATION_REFERENCE] != pmids

    elif isinstance(pmids, Iterable):
        pmids = set(pmids)

        @edge_predicate
        def pmid_exclusion_filter(data: EdgeData) -> bool:
            """Pass for edges with PubMed citations matching one of the contained PubMed identifiers.

            :return: If the edge has a PubMed citation with one of the contained PubMed identifiers
            """
            return has_pubmed(data) and data[CITATION][CITATION_REFERENCE] not in pmids

    else:
        raise TypeError

    return pmid_exclusion_filter


def build_author_inclusion_filter(authors: Strings) -> EdgePredicate:
    """Pass only for edges with author information that matches one of the given authors.

    :param authors: The author or list of authors to filter by
    """
    if isinstance(authors, str):
        @edge_predicate
        def author_filter(data: EdgeData) -> bool:
            """Pass only for edges with citations with an author that matches the contained author.

            :return: If the edge has a citation with an author that matches the the contained author
            """
            return has_authors(data) and authors in data[CITATION][CITATION_AUTHORS]

    elif isinstance(authors, Iterable):
        authors = set(authors)

        @edge_predicate
        def author_filter(data: EdgeData) -> bool:
            """Pass only for edges with citations with an author that matches one or more of the contained authors.

            :return: If the edge has a citation with an author that matches the the contained author
            """
            return has_authors(data) and any(
                author in data[CITATION][CITATION_AUTHORS]
                for author in authors
            )

    else:
        raise TypeError

    return author_filter


def node_has_namespace(node: BaseEntity, namespace: str) -> bool:
    """Pass for nodes that have the given namespace."""
    ns = node.get(NAMESPACE)
    return ns is not None and ns == namespace


def node_has_namespaces(node: BaseEntity, namespaces: Set[str]) -> bool:
    """Pass for nodes that have one of the given namespaces."""
    ns = node.get(NAMESPACE)
    return ns is not None and ns in namespaces


def build_source_namespace_filter(namespaces: Strings) -> EdgePredicate:
    """Pass for edges whose source nodes have the given namespace or one of the given namespaces.

    :param namespaces: The namespace or namespaces to filter by
    """
    if isinstance(namespaces, str):
        def source_namespace_filter(_, u: BaseEntity, __, ___) -> bool:
            return node_has_namespace(u, namespaces)

    elif isinstance(namespaces, Iterable):
        namespaces = set(namespaces)

        def source_namespace_filter(_, u: BaseEntity, __, ___) -> bool:
            return node_has_namespaces(u, namespaces)

    else:
        raise TypeError

    return source_namespace_filter


def build_target_namespace_filter(namespaces: Strings) -> EdgePredicate:
    """Only passes for edges whose target nodes have the given namespace or one of the given namespaces.

    :param namespaces: The namespace or namespaces to filter by
    """
    if isinstance(namespaces, str):
        def target_namespace_filter(_, __, v: BaseEntity, ___) -> bool:
            return node_has_namespace(v, namespaces)

    elif isinstance(namespaces, Iterable):
        namespaces = set(namespaces)

        def target_namespace_filter(_, __, v: BaseEntity, ___) -> bool:
            return node_has_namespaces(v, namespaces)

    else:
        raise TypeError

    return target_namespace_filter
