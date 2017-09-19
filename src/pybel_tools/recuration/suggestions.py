# -*- coding: utf-8 -*-

"""This module contains functions for making namespace suggestions"""

import logging

import requests
from requests.compat import quote_plus

__all__ = [
    'get_user_ols_search_url',
    'get_ols_suggestion',
    'get_ols_search',
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
