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
from pybel.struct.filters.edge_predicates import edge_predicate, has_authors, has_pubmed
from pybel.utils import subdict_matches
from ..utils import check_has_annotation

__all__ = [
    'summarize_edge_filter',
    'build_annotation_value_filter',
    'build_edge_data_filter',
    'build_annotation_dict_all_filter',
    'build_annotation_dict_any_filter',
    'build_relation_filter',
    'build_pmid_inclusion_filter',
    'build_pmid_exclusion_filter',
    'build_author_inclusion_filter',
    'edge_has_activity',
    'edge_has_translocation',
    'edge_has_degradation',
    'edge_has_pathology_causal',
]


def summarize_edge_filter(graph, edge_filters):
    """Prints a summary of the number of edges passing a given set of filters

    :param pybel.BELGraph graph: A BEL graph
    :param edge_filters: A filter or list of filters
    :type edge_filters: types.FunctionType or list[types.FunctionType] or tuple[types.FunctionType]
    """
    passed = count_passed_edge_filter(graph, edge_filters)
    print('{}/{} edges passed {}'.format(
        passed, graph.number_of_edges(),
        ', '.join(edge_filter.__name__ for edge_filter in edge_filters)
    ))


def build_annotation_value_filter(annotation, value):
    """Builds a filter that only passes for edges that contain the given annotation and have the given value(s)
    
    :param str annotation: The annotation to filter on
    :param str or iter[str] value: The value or values for the annotation to filter on
    :return: An edge filter function (graph, node, node, key, data) -> bool
    :rtype: types.FunctionType
    """

    if isinstance(value, str):
        @edge_predicate
        def annotation_value_filter(data):
            """Only passes for edges that contain the given annotation and have the given value
    
            :param dict data: The edge data dictionary
            :return: If the edge has the contained annotation with the contained value
            :rtype: bool
            """
            return check_has_annotation(data, annotation) and data[ANNOTATIONS][annotation] == value

        return annotation_value_filter

    elif isinstance(value, Iterable):
        values = set(value)

        @edge_predicate
        def annotation_values_filter(data):
            """Only passes for edges that contain the given annotation and have one of the given values

            :param dict data: The edge data dictionary
            :return: If the edge has the contained annotation and one of the contained values
            :rtype: bool
            """
            return check_has_annotation(data, annotation) and data[ANNOTATIONS][annotation] in values

        return annotation_values_filter


def build_edge_data_filter(annotations, partial_match=True):
    """Builds a filter that keeps edges whose data dictionaries are superdictionaries to the given dictionary

    :param dict annotations: The annotation query dict to match
    :param bool partial_match: Should the query values be used as partial or exact matches? Defaults to :code:`True`.
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
            check_has_annotation(data, key) and data[ANNOTATIONS][key] == value
            for key, values in annotations.items()
            for value in values
        )

    return annotation_dict_filter


def build_relation_filter(relations):
    """Builds a filter that passes only for edges with the given relation
    
    :param relations: A relation or iterable of relations
    :type relations: str or iter[str]
    :return: An edge filter function (graph, node, node, key, data) -> bool
    :rtype: types.FunctionType
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
    :rtype: types.FunctionType
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
    :rtype: types.FunctionType
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
    :rtype: types.FunctionType
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
                author in data[CITATION][CITATION_AUTHORS]
                for author in authors
            )

        return author_filter


def _edge_has_modifier(d, modifier):
    """Checks if the edge has the given modifier

    :param dict d: The edge data dictionary
    :param str modifier: The modifier to check. One of :data:`pybel.constants.ACTIVITY`,
                        :data:`pybel.constants.DEGRADATION`, or :data:`pybel.constants.TRANSLOCATION`.
    :return: Does either the subject or object have the given modifier
    :rtype: bool
    """
    if SUBJECT in d:
        return MODIFIER in d[SUBJECT] and d[SUBJECT][MODIFIER] == modifier
    elif OBJECT in d:
        return MODIFIER in d[OBJECT] and d[OBJECT][MODIFIER] == modifier
    return False


@edge_predicate
def edge_has_activity(d):
    """Checks if the edge contains an activity in either the subject or object

    :param dict d: The edge data dictionary
    :return: If the edge contains an activity in either the subject or object
    :rtype: bool
    """
    return _edge_has_modifier(d, ACTIVITY)


@edge_predicate
def edge_has_translocation(data):
    """Checks if the edge has a translocation in either the subject or object

    :param dict data: The edge data dictionary
    :return: If the edge has a translocation in either the subject or object
    :rtype: bool
    """
    return _edge_has_modifier(data, TRANSLOCATION)


@edge_predicate
def edge_has_degradation(data):
    """Checks if the edge contains a degradation in either the subject or object

    :param dict data: The edge data dictionary
    :return: If the edge contains a degradation in either the subject or object
    :rtype: bool
    """
    return _edge_has_modifier(data, DEGRADATION)


@edge_predicate
def edge_has_pathology_causal(graph, u, v, k, d):
    """Returns if the subject of this edge is a pathology and participates in a causal relation where the object is
    not a pathology. These relations are generally nonsense.

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: A BEL node
    :param tuple v: A BEL node
    :param int k: The edge key between the given nodes
    :param dict d: The edge data dictionary
    :return: If the subject of this edge is a pathology and it participates in a causal reaction.
    :rtype: bool
    """
    return (
        graph.node[u][FUNCTION] == PATHOLOGY and
        d[RELATION] in CAUSAL_RELATIONS and
        graph.node[v][FUNCTION] not in {PATHOLOGY, BIOPROCESS}
    )
