# -*- coding: utf-8 -*-

"""
Edge Filters
------------

A edge filter is a function that takes five arguments: a :class:`pybel.BELGraph`, a source node tuple, a target node
tuple, a key, and a data dictionary. It returns a boolean representing whether the edge passed the given test.

This module contains a set of default functions for filtering lists of edges and building edge filtering functions.

A general use for an edge filter function is to use the built-in :func:`filter` in code like
:code:`filter(your_edge_filter, graph.edges_iter(keys=True, data=True))`
"""

from collections import Iterable

from pybel import BELGraph
from pybel.constants import *
from pybel.struct.filters import (
    build_annotation_dict_all_filter, build_annotation_dict_any_filter, count_passed_edge_filter,
)
from pybel.struct.filters.edge_predicates import edge_predicate, has_authors, has_pathology_causal, has_pubmed
from pybel.utils import subdict_matches
from .typing import EdgePredicate, EdgePredicates, Strings

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


def summarize_edge_filter(graph: BELGraph, edge_filters: EdgePredicates):
    """Prints a summary of the number of edges passing a given set of filters

    :param graph: A BEL graph
    :param edge_filters: A filter or list of filters
    """
    passed = count_passed_edge_filter(graph, edge_filters)
    print('{}/{} edges passed {}'.format(
        passed, graph.number_of_edges(),
        (
            ', '.join(edge_filter.__name__ for edge_filter in edge_filters) if isinstance(edge_filters, Iterable) else
            edge_filters.__name__
        )
    ))


def build_edge_data_filter(annotations, partial_match=True) -> EdgePredicate:
    """Build a filter that keeps edges whose data dictionaries are super-dictionaries to the given dictionary.

    :param dict annotations: The annotation query dict to match
    :param bool partial_match: Should the query values be used as partial or exact matches? Defaults to :code:`True`.
    """

    @edge_predicate
    def annotation_dict_filter(data):
        """A filter that matches edges with the given dictionary as a sub-dictionary."""
        return subdict_matches(data, annotations, partial_match=partial_match)

    return annotation_dict_filter


def build_pmid_inclusion_filter(pmid: Strings) -> EdgePredicate:
    """Pass for edges with citations whose references are one of the given PubMed identifiers.
    
    :param pmid: A PubMed identifier or list of PubMed identifiers to filter for
    """

    if isinstance(pmid, str):
        @edge_predicate
        def pmid_inclusion_filter(data):
            """Only passes for edges with PubMed citations matching the contained PubMed identifier

            :param dict data: The edge data dictionary
            :return: If the edge has a PubMed citation with the contained PubMed identifier
            :rtype: bool
            """
            return has_pubmed(data) and data[CITATION][CITATION_REFERENCE] == pmid

        return pmid_inclusion_filter

    else:
        pmids = set(pmid)

        @edge_predicate
        def pmid_inclusion_filter(data):
            """Only passes for edges with PubMed citations matching one of the contained PubMed identifiers

            :param dict data: The edge data dictionary
            :return: If the edge has a PubMed citation with one of the contained PubMed identifiers
            :rtype: bool
            """
            return has_pubmed(data) and data[CITATION][CITATION_REFERENCE] in pmids

        return pmid_inclusion_filter


def build_pmid_exclusion_filter(pmid: Strings) -> EdgePredicate:
    """Fail for edges with citations whose references are one of the given PubMed identifiers.

    :param pmid: A PubMed identifier or list of PubMed identifiers to filter against
    """

    if isinstance(pmid, str):

        @edge_predicate
        def pmid_exclusion_filter(data):
            """Fails for edges with PubMed citations matching the contained PubMed identifier

            :param dict data: The edge data dictionary
            :return: If the edge has a PubMed citation with the contained PubMed identifier
            :rtype: bool
            """
            return has_pubmed(data) and data[CITATION][CITATION_REFERENCE] != pmid

        return pmid_exclusion_filter

    else:
        pmids = set(pmid)

        @edge_predicate
        def pmid_exclusion_filter(data):
            """Only passes for edges with PubMed citations matching one of the contained PubMed identifiers

            :param dict data: The edge data dictionary
            :return: If the edge has a PubMed citation with one of the contained PubMed identifiers
            :rtype: bool
            """
            return has_pubmed(data) and data[CITATION][CITATION_REFERENCE] not in pmids

        return pmid_exclusion_filter


def build_author_inclusion_filter(author: Strings) -> EdgePredicate:
    """Only passes for edges with author information that matches one of the given authors
    
    :param author: The author or list of authors to filter by
    """
    if isinstance(author, str):
        @edge_predicate
        def author_filter(data):
            """Only passes for edges with citations with an author that matches the contained author

            :param dict data: The edge data dictionary
            :return: If the edge has a citation with an author that matches the the contained author
            :rtype: bool
            """
            return has_authors(data) and author in data[CITATION][CITATION_AUTHORS]

        return author_filter

    else:
        authors = set(author)

        @edge_predicate
        def author_filter(data):
            """Only passes for edges with citations with an author that matches one or more of the contained authors

            :param dict data: The edge data dictionary
            :return: If the edge has a citation with an author that matches the the contained author
            :rtype: bool
            """
            return has_authors(data) and any(
                a in data[CITATION][CITATION_AUTHORS]
                for a in authors
            )

        return author_filter


def node_has_namespace(graph, node, namespace):
    ns = graph.node[node].get(NAMESPACE)
    return ns is not None and ns == namespace


def node_has_namespaces(graph, node, namespaces):
    namespace = graph.node[node].get(NAMESPACE)
    return namespace is not None and namespace in namespaces


def build_source_namespace_filter(namespaces: Strings) -> EdgePredicate:
    """Only passes for edges whose source nodes have the given namespace or one of the given namespaces

    :param namespaces: The namespace or namespaces to filter by
    """
    if isinstance(namespaces, str):

        def source_namespace_filter(graph, u, v, k):
            return node_has_namespace(graph, u, namespaces)

    elif isinstance(namespaces, Iterable):
        namespaces = set(namespaces)

        def source_namespace_filter(graph, u, v, k):
            return node_has_namespaces(graph, u, namespaces)

    else:
        raise ValueError

    return source_namespace_filter


def build_target_namespace_filter(namespaces: Strings) -> EdgePredicate:
    """Only passes for edges whose target nodes have the given namespace or one of the given namespaces

    :param namespaces: The namespace or namespaces to filter by
    """
    if isinstance(namespaces, str):
        def target_namespace_filter(graph, u, v, k):
            return node_has_namespace(graph, v, namespaces)

    elif isinstance(namespaces, Iterable):
        namespaces = set(namespaces)

        def target_namespace_filter(graph, u, v, k):
            return node_has_namespaces(graph, v, namespaces)

    else:
        raise ValueError

    return target_namespace_filter
