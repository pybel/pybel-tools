# -*- coding: utf-8 -*-

"""This module contains functions for making namespace suggestions"""

import logging

import requests
from requests.compat import quote_plus

__all__ = [
    'get_user_ols_search_url',
    'get_ols_suggestion',
    'get_ols_search',
    'help_suggest_name',
]

log = logging.getLogger(__name__)

OLS_USER_SEARCH_FMT = 'http://www.ebi.ac.uk/ols/search'
OLS_MACHINE_SUGGESTION_FMT = 'http://www.ebi.ac.uk/ols/api/suggest'
OLS_MACHINE_SEARCH_FMT = 'http://www.ebi.ac.uk/ols/api/search'


def get_user_ols_search_url(name):
    """Gets the URL of the page a user should check when they're not sure about an entity's name"""
    return '{}/?q={}'.format(OLS_USER_SEARCH_FMT, quote_plus(name))


def get_ols_suggestion_url(name):
    return '{}/?q={}'.format(OLS_MACHINE_SUGGESTION_FMT, quote_plus(name))


def get_ols_search_url(name):
    return '{}/?q={}'.format(OLS_MACHINE_SEARCH_FMT, quote_plus(name))


def get_ols_suggestion(name, ontology=None):
    """Gets suggestions from the Ontology Lookup Service for which name is best

    :param str name: The name to search
    :param list[str] ontology: Restrict a search to a set of ontologies e.g. ontology=uberon,ma
    :return: The JSON response
    :rtype: dict

    .. seealso:: https://www.ebi.ac.uk/ols/docs/api#_suggest_term
    """
    params = {
        'q': name
    }

    if ontology:
        params['ontology'] = ','.join(ontology)

    response = requests.get(OLS_MACHINE_SUGGESTION_FMT, params=params)
    return response.json()


def get_ols_search(name):
    """Performs a search with the Ontology Lookup Service"""
    res = requests.get(OLS_MACHINE_SEARCH_FMT, params={'q': name})
    return res.json()


def help_suggest_name(namespace, name, metadata_parser, suggestion_cache):
    """Helps populate a suggestion cache for missing names

    :param str namespace: The namespace to search
    :param str name: The putative name in the namespace
    :param metadata_parser: A metadata parser, which contains the namespace dictionary
    :type metadata_parser: pybel.parser.parse_metadata.MetadataParser
    :param dict suggestion_cache: A defaultdict of lists
    :return:
    :rtype: dict
    """
    from fuzzywuzzy import process, fuzz

    if (namespace, name) in suggestion_cache:
        return suggestion_cache[namespace, name]

    if namespace not in metadata_parser.namespace_dict:
        raise ValueError('Namespace not cached: {}'.format(namespace))

    terms = set(metadata_parser.namespace_dict[namespace])

    for putative, _ in process.extract(name, terms, scorer=fuzz.partial_token_sort_ratio, limit=5):
        suggestion_cache[namespace, name].append(putative)

    return suggestion_cache[namespace, name]
