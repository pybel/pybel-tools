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

from pybel.constants import *
from pybel.struct.filters import count_passed_edge_filter
from pybel.struct.filters.edge_predicates import (
    edge_has_annotation, edge_predicate, has_authors, has_pubmed, is_causal_relation,
)
from pybel.struct.filters.node_predicates import is_pathology
from pybel.utils import subdict_matches

__all__ = [
    'summarize_edge_filter',
    'build_edge_data_filter',
    'build_annotation_dict_all_filter',
    'build_annotation_dict_any_filter',
    'build_relation_filter',
    'build_pmid_inclusion_filter',
    'build_pmid_exclusion_filter',
    'build_author_inclusion_filter',
    'has_pathology_causal',
    'build_source_namespace_filter',
    'build_target_namespace_filter',
]


def summarize_edge_filter(graph, edge_filters):
    """Prints a summary of the number of edges passing a given set of filters

    :param pybel.BELGraph graph: A BEL graph
    :param edge_filters: A filter or list of filters
    :type edge_filters: (pybel.BELGraph, tuple, tuple, int) -> bool or iter[(pybel.BELGraph, tuple, tuple, int) -> bool]
    """
    passed = count_passed_edge_filter(graph, edge_filters)
    print('{}/{} edges passed {}'.format(
        passed, graph.number_of_edges(),
        (
            ', '.join(edge_filter.__name__ for edge_filter in edge_filters) if isinstance(edge_filters, Iterable) else
            edge_filters.__name__
        )
    ))


def build_edge_data_filter(annotations, partial_match=True):
    """Builds a filter that keeps edges whose data dictionaries are superdictionaries to the given dictionary

    :param dict annotations: The annotation query dict to match
    :param bool partial_match: Should the query values be used as partial or exact matches? Defaults to :code:`True`.
    :rtype: (pybel.BELGraph, tuple, tuple, int) -> bool
    """

    @edge_predicate
    def annotation_dict_filter(data):
        """A filter that matches edges with the given dictionary as a subdictionary"""
        return subdict_matches(data, annotations, partial_match=partial_match)

    return annotation_dict_filter


def build_annotation_dict_all_filter(annotations, partial_match=True):
    """Builds a filter that keeps edges whose data dictionaries's annotations entry are superdictionaries to the given
    dictionary

    :param dict annotations: The annotation query dict to match
    :param bool partial_match: Should the query values be used as partial or exact matches? Defaults to :code:`True`.
    :rtype: (pybel.BELGraph, tuple, tuple, int) -> bool
    """

    @edge_predicate
    def annotation_dict_filter(data):
        """A filter that matches edges with the given dictionary as a subdictionary"""
        return subdict_matches(data[ANNOTATIONS], annotations, partial_match=partial_match)

    return annotation_dict_filter


def build_annotation_dict_any_filter(annotations):
    """Builds a filter that keeps edges whose data dictionaries's annotations entry contain any match to
    the target dictionary

    :param dict annotations: The annotation query dict to match
    """

    @edge_predicate
    def annotation_dict_filter(data):
        """A filter that matches edges with the given dictionary as a subdictionary

        :param dict data: A PyBEL edge data dictionary
        :rtype: bool
        """
        return any(
            edge_has_annotation(data, key) and data[ANNOTATIONS][key] == value
            for key, values in annotations.items()
            for value in values
        )

    return annotation_dict_filter


def build_relation_filter(relations):
    """Builds a filter that passes only for edges with the given relation
    
    :param relations: A relation or iterable of relations
    :type relations: str or iter[str]
    :return: An edge filter function (graph, node, node, key, data) -> bool
    :rtype: (pybel.BELGraph, tuple, tuple, int) -> bool
    """

    if isinstance(relations, str):
        @edge_predicate
        def relation_filter(data):
            """Only passes for edges with the contained relation

            :param dict data: A PyBEL edge data dictionary
            :return: If the edge has the contained relation
            :rtype: bool
            """
            return data[RELATION] == relations

        return relation_filter

    elif isinstance(relations, Iterable):
        relation_set = set(relations)

        @edge_predicate
        def relation_filter(data):
            """Only passes for edges with one of the contained relations

            :param dict data: A PyBEL edge data dictionary
            :return: If the edge has one of the contained relations
            :rtype: bool
            """
            return data[RELATION] in relation_set

        return relation_filter

    else:
        raise ValueError('Invalid type for argument: {}'.format(relations))


def build_pmid_inclusion_filter(pmid):
    """Only passes for edges with citations whose references are one of the given PubMed identifiers
    
    :param pmid: A PubMed identifier or list of PubMed identifiers to filter for
    :type pmid: str or iter[str]
    :return: An edge filter function (graph, node, node, key, data) -> bool
    :rtype: (pybel.BELGraph, tuple, tuple, int) -> bool
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


def build_pmid_exclusion_filter(pmid):
    """Fails for edges with citations whose references are one of the given PubMed identifiers

    :param pmid: A PubMed identifier or list of PubMed identifiers to filter against
    :type pmid: str or iter[str]
    :return: An edge filter function (graph, node, node, key, data) -> bool
    :rtype: (pybel.BELGraph, tuple, tuple, int) -> bool
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


def build_author_inclusion_filter(author):
    """Only passes for edges with author information that matches one of the given authors
    
    :param author: The author or list of authors to filter by
    :type author: str or iter[str]
    :return: An edge filter
    :rtype: (pybel.BELGraph, tuple, tuple, int) -> bool
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


def has_pathology_causal(graph, u, v, k):
    """Returns if the subject of this edge is a pathology and participates in a causal relation where the object is
    not a pathology. These relations are generally nonsense.

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: A BEL node
    :param tuple v: A BEL node
    :param int k: The edge key between the given nodes
    :return: If the subject of this edge is a pathology and it participates in a causal reaction.
    :rtype: bool
    """
    return (
            is_pathology(graph, u) and
            is_causal_relation(graph, u, v, k) and
            graph.node[v][FUNCTION] not in {PATHOLOGY, BIOPROCESS}
    )


def node_has_namespace(graph, node, namespace):
    ns = graph.node[node].get(NAMESPACE)
    return ns is not None and ns == namespace


def node_has_namespaces(graph, node, namespaces):
    namespace = graph.node[node].get(NAMESPACE)
    return namespace is not None and namespace in namespaces


def build_source_namespace_filter(namespaces):
    """Only passes for edges whose source nodes have the given namespace or one of the given namespaces

    :param namespaces: The namespace or namespaces to filter by
    :type namespaces: str or iter[str]
    :rtype: (pybel.BELGraph, tuple, tuple, int) -> bool
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


def build_target_namespace_filter(namespaces):
    """Only passes for edges whose target nodes have the given namespace or one of the given namespaces

    :param namespaces: The namespace or namespaces to filter by
    :type namespaces: str or iter[str]
    :rtype: (pybel.BELGraph, tuple, tuple, int) -> bool
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
