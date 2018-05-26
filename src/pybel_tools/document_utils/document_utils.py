# -*- coding: utf-8 -*-

"""Utilities to merge multiple BEL documents on the same topic"""

import logging
from xml.etree import ElementTree

import requests

from pybel.resources.constants import *
from pybel.resources.document import make_knowledge_header
from pybel_tools.constants import abstract_url_fmt, title_url_fmt

__all__ = [
    'write_boilerplate',
]

log = logging.getLogger(__name__)


def make_pubmed_abstract_group(pmids):
    """Builds a skeleton for the citations' statements
    
    :param pmids: A list of PubMed identifiers
    :type pmids: iter[str] or iter[int]
    :return: An iterator over the lines of the citation section
    :rtype: iter[str]
    """
    for pmid in set(pmids):
        yield ''

        res = requests.get(title_url_fmt.format(pmid))
        title = res.content.decode('utf-8').strip()

        yield citation_format.format(title, pmid)

        res = requests.get(abstract_url_fmt.format(pmid))
        abstract = res.content.decode('utf-8').strip()

        yield evidence_format.format(abstract)
        yield '\nUNSET Evidence\nUNSET Citation'


def _sanitize(s):
    if s is None:
        return None
    return s.strip().replace('\n', '')


#: Allows for querying the Entrez Gene Summary utility by formatting with an entrez id or list of comma seperated ids
PUBMED_GENE_QUERY_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={}'


def get_entrez_gene_data(entrez_ids):
    """Gets gene info from Entrez"""
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


def make_pubmed_gene_group(entrez_ids):
    """Builds a skeleton for gene summaries

    :param list[str] entrez_ids: A list of entrez id's to query the pubmed service 
    :return: An iterator over statement lines for NCBI entrez gene summaries
    :rtype: iter[str]
    """
    url = PUBMED_GENE_QUERY_URL.format(','.join(str(x).strip() for x in entrez_ids))
    response = requests.get(url)
    tree = ElementTree.fromstring(response.content)

    for x in tree.findall('./DocumentSummarySet/DocumentSummary'):
        yield '\n# {}'.format(x.find('Description').text)
        yield 'SET Citation = {{"Other", "PubMed Gene", "{}"}}'.format(x.attrib['uid'])
        yield 'SET Evidence = "{}"'.format(x.find('Summary').text.strip().replace('\n', ''))
        yield '\nUNSET Evidence\nUNSET Citation'


def write_boilerplate(name, version=None, description=None, authors=None, contact=None, copyright=None, licenses=None,
                      disclaimer=None, namespace_url=None, namespace_patterns=None, annotation_url=None,
                      annotation_patterns=None, annotation_list=None, pmids=None, entrez_ids=None, file=None):
    """Writes a boilerplate BEL document, with standard document metadata, definitions.

    :param str name: The unique name for this BEL document
    :param str contact: The email address of the maintainer
    :param str description: A description of the contents of this document
    :param str authors: The authors of this document
    :param str version: The version. Defaults to current date in format ``YYYYMMDD``.
    :param str copyright: Copyright information about this document
    :param str licenses: The license applied to this document
    :param str disclaimer: The disclaimer for this document
    :param dict[str,str] namespace_url: an optional dictionary of {str name: str URL} of namespaces
    :param dict[str,str] namespace_patterns: An optional dictionary of {str name: str regex} namespaces
    :param dict[str,str] annotation_url: An optional dictionary of {str name: str URL} of annotations
    :param dict[str,str] annotation_patterns: An optional dictionary of {str name: str regex} of regex annotations
    :param dict[str,set[str]] annotation_list: An optional dictionary of {str name: set of names} of list annotations
    :param iter[str] or iter[int] pmids: A list of PubMed identifiers to auto-populate with citation and abstract
    :param iter[str] or iter[int] entrez_ids: A list of Entrez identifiers to autopopulate the gene summary as evidence
    :param file file: A writable file or file-like. If None, defaults to :data:`sys.stdout`
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
