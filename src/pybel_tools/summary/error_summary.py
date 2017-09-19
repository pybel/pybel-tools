# -*- coding: utf-8 -*-

"""This module contains functions that provide summaries of the errors encountered while parsing a BEL script"""

from collections import Counter, defaultdict

from pybel.constants import ANNOTATIONS
from pybel.parser.parse_exceptions import *
from .node_summary import get_namespaces, get_names_by_namespace
from ..utils import check_has_annotation

__all__ = [
    'count_error_types',
    'count_naked_names',
    'get_naked_names',
    'get_incorrect_names_by_namespace',
    'get_incorrect_names',
    'get_undefined_namespaces',
    'get_undefined_namespace_names',
    'calculate_incorrect_name_dict',
    'calculate_error_by_annotation',
    'group_errors',
    'get_names_including_errors',
    'get_names_including_errors_by_namespace',
    'get_undefined_annotations',
    'get_namespaces_with_incorrect_names',
]


def count_error_types(graph):
    """Counts the occurrence of each type of error in a graph

    :param pybel.BELGraph graph: A BEL graph
    :return: A Counter of {error type: frequency}
    :rtype: collections.Counter
    """
    return Counter(e.__class__.__name__ for _, _, e, _ in graph.warnings)


def _naked_names_iter(graph):
    """Iterates over naked name warnings froma  graph

    :param pybel.BELGraph graph: A BEL graph
    :rtype: iter[NakedNameWarning]
    """
    for _, _, e, _ in graph.warnings:
        if isinstance(e, NakedNameWarning):
            yield e.name


def count_naked_names(graph):
    """Counts the frequency of each naked name (names without namespaces)

    :param pybel.BELGraph graph: A BEL graph
    :return: A Counter from {name: frequency}
    :rtype: collections.Counter
    """
    return Counter(_naked_names_iter(graph))


def get_naked_names(graph):
    """Gets the set of naked names in the graph

    :param pybel.BELGraph graph: A BEL graph
    :rtype: set[str]
    """
    return set(_naked_names_iter(graph))


def get_namespaces_with_incorrect_names(graph):
    """Returns the set of all namespaces with incorrect names in the graph"""
    return {
        e.namespace
        for _, _, e, _ in graph.warnings
        if isinstance(e, (MissingNamespaceNameWarning, MissingNamespaceRegexWarning))
    }


def get_incorrect_names_by_namespace(graph, namespace):
    """Returns the set of all incorrect names from the given namespace in the graph

    :param pybel.BELGraph graph: A BEL graph
    :param str namespace: The namespace to filter by
    :return: The set of all incorrect names from the given namespace in the graph
    :rtype: set[str]
    """
    return {
        e.name
        for _, _, e, _ in graph.warnings
        if isinstance(e, (MissingNamespaceNameWarning, MissingNamespaceRegexWarning)) and e.namespace == namespace
    }


def get_incorrect_names(graph):
    """Returns the dict of the sets of all incorrect names from the given namespace in the graph

    :param pybel.BELGraph graph: A BEL graph
    :return: The set of all incorrect names from the given namespace in the graph
    :rtype: dict[str,set[str]]
    """
    return {
        namespace: get_incorrect_names_by_namespace(graph, namespace)
        for namespace in get_namespaces(graph)
    }


def get_undefined_namespaces(graph):
    """Gets all namespaces that aren't actually defined
    
    :param pybel.BELGraph graph: A BEL graph
    :return: The set of all undefined namespaces
    :rtype: set[str]
    """
    return {
        e.namespace
        for _, _, e, _ in graph.warnings
        if isinstance(e, UndefinedNamespaceWarning)
    }


def get_undefined_namespace_names(graph, namespace):
    """Gets the names from a namespace that wasn't actually defined
    
    :param pybel.BELGraph graph: A BEL graph
    :param str namespace: The namespace to filter by
    :return: The set of all names from the undefined namespace
    :rtype: set[str]
    """
    return {
        e.name
        for _, _, e, _ in graph.warnings
        if isinstance(e, UndefinedNamespaceWarning) and e.namespace == namespace
    }


def get_undefined_annotations(graph):
    """Gets all annotations that aren't actually defined
    
    :param pybel.BELGraph graph: A BEL graph
    :return: The set of all undefined annotations
    :rtype: set[str]
    """
    return {
        e.annotation
        for _, _, e, _ in graph.warnings
        if isinstance(e, UndefinedAnnotationWarning)
    }


# FIXME need to change underlying definition and usage of this exception
def get_undefined_annotation_values(graph, annotation):
    """Gets the values from an annotation that wasn't actually defined

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotaton to filter by
    :return: The set of all values from the undefined annotation
    :rtype: set[str]
    """
    raise NotImplementedError
    # return {e.value for _, _, e, _ in graph.warnings if isinstance(e, UndefinedAnnotationWarning) and e.annotation == annotation}


def calculate_incorrect_name_dict(graph):
    """Groups all of the incorrect identifiers in a dict of {namespace: list of erroneous names}

    :param pybel.BELGraph graph: A BEL graph
    :return: A dictionary of {namespace: list of erroneous names}
    :rtype: dict[str, str]
    """
    missing = defaultdict(list)

    for line_number, line, e, ctx in graph.warnings:
        if not isinstance(e, (MissingNamespaceNameWarning, MissingNamespaceRegexWarning)):
            continue
        missing[e.namespace].append(e.name)

    return dict(missing)


def calculate_error_by_annotation(graph, annotation):
    """Groups the graph by a given annotation and builds lists of errors for each

    :param pybel.BELGraph graph: A BEL graph
    :param annotation: The annotation to group errors by
    :type annotation: str
    :return: A dictionary of {annotation value: list of errors}
    :rtype: dict[str, list[str]]
    """
    results = defaultdict(list)

    for line_number, line, e, context in graph.warnings:
        if not context or not check_has_annotation(context, annotation):
            continue

        values = context[ANNOTATIONS][annotation]

        if isinstance(values, str):
            results[values].append(e.__class__.__name__)
        elif isinstance(values, (set, tuple, list)):
            for value in values:
                results[value].append(e.__class__.__name__)

    return dict(results)


def group_errors(graph):
    """Groups the errors together for analysis of the most frequent error

    :param pybel.BELGraph graph: A BEL graph
    :return: A dictionary of {error string: list of line numbers}
    :rtype: dict[str, list[int]]
    """
    warning_summary = defaultdict(list)

    for ln, _, e, _ in graph.warnings:
        warning_summary[str(e)].append(ln)

    return dict(warning_summary)


def get_names_including_errors_by_namespace(graph, namespace):
    """Takes the names from the graph in a given namespace and the erroneous names from the same namespace and returns
    them together as a unioned set

    :param pybel.BELGraph graph: A BEL graph
    :param str namespace: The namespace to filter by
    :return: The set of all correct and incorrect names from the given namespace in the graph
    :rtype: set[str]
    """
    return get_names_by_namespace(graph, namespace) | get_incorrect_names_by_namespace(graph, namespace)


def get_names_including_errors(graph):
    """Takes the names from the graph in a given namespace and the erroneous names from the same namespace and returns
    them together as a unioned set

    :param pybel.BELGraph graph: A BEL graph
    :return: The dict of the sets of all correct and incorrect names from the given namespace in the graph
    :rtype: dict[str,set[str]]
    """
    return {
        namespace: get_names_including_errors_by_namespace(graph, namespace)
        for namespace in get_namespaces(graph)
    }
