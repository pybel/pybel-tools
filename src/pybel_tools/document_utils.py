# -*- coding: utf-8 -*-

"""Utilities to merge multiple BEL documents on the same topic."""

import logging
from typing import Iterable, Mapping, Optional, Set, TextIO, Union
from xml.etree import ElementTree

import requests

from bel_resources import make_knowledge_header

__all__ = [
    'write_boilerplate',
]

log = logging.getLogger(__name__)

abstract_url_fmt = "http://togows.dbcls.jp/entry/ncbi-pubmed/{}/abstract"
title_url_fmt = "http://togows.dbcls.jp/entry/ncbi-pubmed/{}/title"
#: SO gives short citation information
so_url_fmt = "http://togows.dbcls.jp/entry/ncbi-pubmed/{}/so"


def make_pubmed_abstract_group(pmids: Iterable[Union[str, int]]) -> Iterable[str]:
    """Build a skeleton for the citations' statements.

    :param pmids: A list of PubMed identifiers
    :return: An iterator over the lines of the citation section
    """
    for pmid in set(pmids):
        yield ''

        res = requests.get(title_url_fmt.format(pmid))
        title = res.content.decode('utf-8').strip()

        yield f'SET Citation = {{"{title}", "{pmid}"}}'

        res = requests.get(abstract_url_fmt.format(pmid))
        abstract = res.content.decode('utf-8').strip()

        yield f'SET Evidence = "{abstract}"'
        yield '\nUNSET Evidence\nUNSET Citation'


def _sanitize(s):
    if s is not None:
        return s.strip().replace('\n', '')


#: Allows for querying the Entrez Gene Summary utility by formatting with an entrez id or list of comma seperated ids
PUBMED_GENE_QUERY_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={}'


def get_entrez_gene_data(entrez_ids: Iterable[Union[str, int]]):
    """Get gene info from Entrez."""
    url = PUBMED_GENE_QUERY_URL.format(','.join(str(x).strip() for x in entrez_ids))
    response = requests.get(url)
    tree = ElementTree.fromstring(response.content)

    return {
        element.attrib['uid']: {
            'summary': _sanitize(element.find('Summary').text),
            'description': element.find('Description').text
        }
        for element in tree.findall('./DocumentSummarySet/DocumentSummary')
    }


def make_pubmed_gene_group(entrez_ids: Iterable[Union[str, int]]) -> Iterable[str]:
    """Build a skeleton for gene summaries.

    :param entrez_ids: A list of Entrez Gene identifiers to query the PubMed service
    :return: An iterator over statement lines for NCBI Entrez Gene summaries
    """
    url = PUBMED_GENE_QUERY_URL.format(','.join(str(x).strip() for x in entrez_ids))
    response = requests.get(url)
    tree = ElementTree.fromstring(response.content)

    for x in tree.findall('./DocumentSummarySet/DocumentSummary'):
        yield '\n# {}'.format(x.find('Description').text)
        yield 'SET Citation = {{"Other", "PubMed Gene", "{}"}}'.format(x.attrib['uid'])
        yield 'SET Evidence = "{}"'.format(x.find('Summary').text.strip().replace('\n', ''))
        yield '\nUNSET Evidence\nUNSET Citation'


def write_boilerplate(
        name: str,
        version: Optional[str] = None,
        description: Optional[str] = None,
        authors: Optional[str] = None,
        contact: Optional[str] = None,
        copyright: Optional[str] = None,
        licenses: Optional[str] = None,
        disclaimer: Optional[str] = None,
        namespace_url: Optional[Mapping[str, str]] = None,
        namespace_patterns: Optional[Mapping[str, str]] = None,
        annotation_url: Optional[Mapping[str, str]] = None,
        annotation_patterns: Optional[Mapping[str, str]] = None,
        annotation_list: Optional[Mapping[str, Set[str]]] = None,
        pmids: Optional[Iterable[Union[str, int]]] = None,
        entrez_ids: Optional[Iterable[Union[str, int]]] = None,
        file: Optional[TextIO] = None,
) -> None:
    """Write a boilerplate BEL document, with standard document metadata, definitions.

    :param name: The unique name for this BEL document
    :param contact: The email address of the maintainer
    :param description: A description of the contents of this document
    :param authors: The authors of this document
    :param version: The version. Defaults to current date in format ``YYYYMMDD``.
    :param copyright: Copyright information about this document
    :param licenses: The license applied to this document
    :param disclaimer: The disclaimer for this document
    :param namespace_url: an optional dictionary of {str name: str URL} of namespaces
    :param namespace_patterns: An optional dictionary of {str name: str regex} namespaces
    :param annotation_url: An optional dictionary of {str name: str URL} of annotations
    :param annotation_patterns: An optional dictionary of {str name: str regex} of regex annotations
    :param annotation_list: An optional dictionary of {str name: set of names} of list annotations
    :param pmids: A list of PubMed identifiers to auto-populate with citation and abstract
    :param entrez_ids: A list of Entrez identifiers to autopopulate the gene summary as evidence
    :param file: A writable file or file-like. If None, defaults to :data:`sys.stdout`
    """
    lines = make_knowledge_header(
        name=name,
        version=version or '1.0.0',
        description=description,
        authors=authors,
        contact=contact,
        copyright=copyright,
        licenses=licenses,
        disclaimer=disclaimer,
        namespace_url=namespace_url,
        namespace_patterns=namespace_patterns,
        annotation_url=annotation_url,
        annotation_patterns=annotation_patterns,
        annotation_list=annotation_list,
    )

    for line in lines:
        print(line, file=file)

    if pmids is not None:
        for line in make_pubmed_abstract_group(pmids):
            print(line, file=file)

    if entrez_ids is not None:
        for line in make_pubmed_gene_group(entrez_ids):
            print(line, file=file)
