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

from pybel.constants import *
from pybel.struct.filters import count_passed_edge_filter
from pybel.utils import subdict_matches
from ..constants import PUBMED
from ..utils import check_has_annotation

__all__ = [
    'edge_is_causal',
    'edge_has_author_annotation',
    'build_inverse_filter',
    'build_annotation_value_filter',
    'build_edge_data_filter',
    'build_pmid_inclusion_filter',
    'build_pmid_exclusion_filter',
    'build_author_inclusion_filter',
    'build_relation_filter',
    'summarize_edge_filter',
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


def edge_is_causal(graph, u, v, k, d):
    """Only passes on causal edges, belonging to the set :data:`pybel.constants.CAUSAL_RELATIONS`

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: A BEL node
    :param tuple v: A BEL node
    :param int k: The edge key between the given nodes
    :param dict d: The edge data dictionary
    :return: If the edge is a causal edge
    :rtype: bool
    """
    return graph.edge[u][v][k][RELATION] in CAUSAL_RELATIONS


def edge_has_author_annotation(graph, u, v, k, d):
    """Passes for edges that have citations with authors

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: A BEL node
    :param tuple v: A BEL node
    :param int k: The edge key between the given nodes
    :param dict d: The edge data dictionary
    :return: Does the edge's citation data dictionary have authors included?
    :rtype: bool
    """
    return CITATION in graph.edge[u][v][k] and CITATION_AUTHORS in graph.edge[u][v][k][CITATION]


def edge_has_pubmed_citation(graph, u, v, k, d):
    """Passes for edges that have PubMed citations
    
    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: A BEL node
    :param tuple v: A BEL node
    :param int k: The edge key between the given nodes
    :param dict d: The edge data dictionary
    :return: Is the edge's citation from :data:`PUBMED`?
    :rtype: bool
    """
    return CITATION in graph.edge[u][v][k] and PUBMED == graph.edge[u][v][k][CITATION][CITATION_TYPE]


def build_inverse_filter(edge_filter):
    """Builds a filter that is the inverse of the given filter
    
    :param edge_filter: An edge filter function (graph, node, node, key, data) -> bool
    :type edge_filter: types.FunctionType
    :return: An edge filter function (graph, node, node, key, data) -> bool
    :rtype: types.FunctionType
    """

    def inverse_filter(graph, u, v, k, d):
        return not edge_filter(graph, u, v, k, d)

    return inverse_filter


def build_annotation_value_filter(annotation, value):
    """Builds a filter that only passes for edges that contain the given annotation and have the given value(s)
    
    :param str annotation: The annotation to filter on
    :param str or iter[str] value: The value or values for the annotation to filter on
    :return: An edge filter function (graph, node, node, key, data) -> bool
    :rtype: types.FunctionType
    """

    if isinstance(value, str):
        def annotation_value_filter(graph, u, v, k, d):
            """Only passes for edges that contain the given annotation and have the given value
    
            :param pybel.BELGraph graph: A BEL Graph
            :param tuple u: A BEL node
            :param tuple v: A BEL node
            :param int k: The edge key between the given nodes
            :param dict d: The edge data dictionary
            :return: If the edge has the contained annotation with the contained value
            :rtype: bool
            """
            if not check_has_annotation(graph.edge[u][v][k], annotation):
                return False

            return graph.edge[u][v][k][ANNOTATIONS][annotation] == value

        return annotation_value_filter

    elif isinstance(value, (list, set, tuple)):
        values = set(value)

        def annotation_values_filter(graph, u, v, k, d):
            """Only passes for edges that contain the given annotation and have one of the given values

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple u: A BEL node
            :param tuple v: A BEL node
            :param int k: The edge key between the given nodes
            :param dict d: The edge data dictionary
            :return: If the edge has the contained annotation and one of the contained values
            :rtype: bool
            """
            if not check_has_annotation(graph.edge[u][v][k], annotation):
                return False

            return graph.edge[u][v][k][ANNOTATIONS][annotation] in values

        return annotation_values_filter


def build_edge_data_filter(annotations, partial_match=True):
    """Builds a filter that keeps edges whose data dictionaries are superdictionaries to the given dictionary

    :param dict annotations: The annotation query dict to match
    :param bool partial_match: Should the query values be used as partial or exact matches? Defaults to :code:`True`.
    """

    def annotation_dict_filter(graph, u, v, k, d):
        """A filter that matches edges with the given dictionary as a subdictionary"""
        return subdict_matches(d, annotations, partial_match=partial_match)

    return annotation_dict_filter


def build_annotation_dict_all_filter(annotations, partial_match=True):
    """Builds a filter that keeps edges whose data dictionaries's annotations entry are superdictionaries to the given
    dictionary

    :param dict annotations: The annotation query dict to match
    :param bool partial_match: Should the query values be used as partial or exact matches? Defaults to :code:`True`.
    """

    def annotation_dict_filter(graph, u, v, k, d):
        """A filter that matches edges with the given dictionary as a subdictionary"""
        return subdict_matches(d[ANNOTATIONS], annotations, partial_match=partial_match)

    return annotation_dict_filter


def build_annotation_dict_any_filter(annotations):
    """Builds a filter that keeps edges whose data dictionaries's annotations entry contain any match to
    the target dictionary

    :param dict annotations: The annotation query dict to match
    """

    def annotation_dict_filter(graph, u, v, k, d):
        """A filter that matches edges with the given dictionary as a subdictionary"""
        return any(
            ANNOTATIONS in d and key in d[ANNOTATIONS] and d[ANNOTATIONS][key] == value
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
        def relation_filter(graph, u, v, k, d):
            """Only passes for edges with the contained relation

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple u: A BEL node
            :param tuple v: A BEL node
            :param int k: The edge key between the given nodes
            :param dict d: The edge data dictionary
            :return: If the edge has the contained relation
            :rtype: bool
            """
            return graph.edge[u][v][k][RELATION] == relations

        return relation_filter

    elif isinstance(relations, (list, tuple, set)):
        relation_set = set(relations)

        def relation_filter(graph, u, v, k, d):
            """Only passes for edges with one of the contained relations

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple u: A BEL node
            :param tuple v: A BEL node
            :param int k: The edge key between the given nodes
            :param dict d: The edge data dictionary
            :return: If the edge has one of the contained relations
            :rtype: bool
            """
            return graph.edge[u][v][k][RELATION] in relation_set

        return relation_filter
    else:
        raise ValueError('Invalid type for argument: {}'.format(relations))


def build_pmid_inclusion_filter(pmids):
    """Only passes for edges with citations whose references are one of the given PubMed identifiers
    
    :param pmids: A PubMed identifier or list of PubMed identifiers to filter for
    :type pmids: str or iter[str]
    :return: An edge filter function (graph, node, node, key, data) -> bool
    :rtype: types.FunctionType
    """

    if isinstance(pmids, str):
        def pmid_inclusion_filter(graph, u, v, k, d):
            """Only passes for edges with PubMed citations matching the contained PubMed identifier

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple u: A BEL node
            :param tuple v: A BEL node
            :param int k: The edge key between the given nodes
            :param dict d: The edge data dictionary
            :return: If the edge has a PubMed citation with the contained PubMed identifier
            :rtype: bool
            """
            if not edge_has_pubmed_citation(graph, u, v, k, d):
                return False

            return d[CITATION][CITATION_REFERENCE] == pmids

        return pmid_inclusion_filter

    else:
        pmids = set(pmids)

        def pmid_inclusion_filter(graph, u, v, k, d):
            """Only passes for edges with PubMed citations matching one of the contained PubMed identifiers

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple u: A BEL node
            :param tuple v: A BEL node
            :param int k: The edge key between the given nodes
            :param dict d: The edge data dictionary
            :return: If the edge has a PubMed citation with one of the contained PubMed identifiers
            :rtype: bool
            """
            if not edge_has_pubmed_citation(graph, u, v, k, d):
                return False

            return d[CITATION][CITATION_REFERENCE] in pmids

        return pmid_inclusion_filter


def build_pmid_exclusion_filter(pmids):
    """Fails for edges with citations whose references are one of the given PubMed identifiers

    :param pmids: A PubMed identifier or list of PubMed identifiers to filter against
    :type pmids: str or iter[str]
    :return: An edge filter function (graph, node, node, key, data) -> bool
    :rtype: types.FunctionType
    """

    if isinstance(pmids, str):
        def pmid_exclusion_filter(graph, u, v, k, d):
            """Fails for edges with PubMed citations matching the contained PubMed identifier

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple u: A BEL node
            :param tuple v: A BEL node
            :param int k: The edge key between the given nodes
            :param dict d: The edge data dictionary
            :return: If the edge has a PubMed citation with the contained PubMed identifier
            :rtype: bool
            """
            if not edge_has_pubmed_citation(graph, u, v, k, d):
                return False

            return d[CITATION][CITATION_REFERENCE] != pmids

        return pmid_exclusion_filter

    else:
        pmids = set(pmids)

        def pmid_exclusion_filter(graph, u, v, k, d):
            """Only passes for edges with PubMed citations matching one of the contained PubMed identifiers

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple u: A BEL node
            :param tuple v: A BEL node
            :param int k: The edge key between the given nodes
            :param dict d: The edge data dictionary
            :return: If the edge has a PubMed citation with one of the contained PubMed identifiers
            :rtype: bool
            """
            if not edge_has_pubmed_citation(graph, u, v, k, d):
                return False
            return d[CITATION][CITATION_REFERENCE] not in pmids

        return pmid_exclusion_filter


def build_author_inclusion_filter(authors):
    """Only passes for edges with author information that matches one of the given authors
    
    :param authors: The author or list of authors to filter by
    :type authors: str or iter[str]
    :return: An edge filter
    :rtype: types.FunctionType
    """
    if isinstance(authors, str):
        def author_filter(graph, u, v, k, d):
            """Only passes for edges with citations with an author that matches the contained author

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple u: A BEL node
            :param tuple v: A BEL node
            :param int k: The edge key between the given nodes
            :param dict d: The edge data dictionary
            :return: If the edge has a citation with an author that matches the the contained author
            :rtype: bool
            """
            if not edge_has_author_annotation(graph, u, v, k, d):
                return False

            return authors in d[CITATION][CITATION_AUTHORS]

        return author_filter

    else:
        author_set = set(authors)

        def author_filter(graph, u, v, k, d):
            """Only passes for edges with citations with an author that matches one or more of the contained authors

            :param pybel.BELGraph graph: A BEL Graph
            :param tuple u: A BEL node
            :param tuple v: A BEL node
            :param int k: The edge key between the given nodes
            :param dict d: The edge data dictionary
            :return: If the edge has a citation with an author that matches the the contained author
            :rtype: bool
            """
            if not edge_has_author_annotation(graph, u, v, k, d):
                return False

            return any(
                author
                in d[CITATION][CITATION_AUTHORS]
                for author in author_set
            )

        return author_filter


def _edge_has_modifier(graph, u, v, k, d, modifier):
    """Checks if the edge has the given modifier

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: A BEL node
    :param tuple v: A BEL node
    :param int k: The edge key between the given nodes
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


def edge_has_activity(graph, u, v, k, d):
    """Checks if the edge contains an activity in either the subject or object

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: A BEL node
    :param tuple v: A BEL node
    :param int k: The edge key between the given nodes
    :param dict d: The edge data dictionary
    :return: If the edge contains an activity in either the subject or object
    :rtype: bool
    """
    return _edge_has_modifier(graph, u, v, k, d, ACTIVITY)


def edge_has_translocation(graph, u, v, k, d):
    """Checks if the edge has a translocation in either the subject or object

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: A BEL node
    :param tuple v: A BEL node
    :param int k: The edge key between the given nodes
    :param dict d: The edge data dictionary
    :return: If the edge has a translocation in either the subject or object
    :rtype: bool
    """
    return _edge_has_modifier(graph, u, v, k, d, TRANSLOCATION)


def edge_has_degradation(graph, u, v, k, d):
    """Checks if the edge contains a degradation in either the subject or object

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: A BEL node
    :param tuple v: A BEL node
    :param int k: The edge key between the given nodes
    :param dict d: The edge data dictionary
    :return: If the edge contains a degradation in either the subject or object
    :rtype: bool
    """
    return _edge_has_modifier(graph, u, v, k, d, DEGRADATION)


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
