# -*- coding: utf-8 -*-

"""Utilities to merge multiple BEL documents on the same topic"""

from __future__ import print_function

import logging
import os
import sys
import time
from itertools import islice
from xml.etree import ElementTree

import requests

from .constants import title_url_fmt, citation_format, abstract_url_fmt, evidence_format
from .resources import default_namespaces, default_namespace_patterns, default_annotations
from .utils import get_version

__all__ = [
    'merge',
    'make_document_metadata',
    'make_document_namespaces',
    'make_document_annotations',
    'make_pubmed_abstract_group',
    'get_entrez_gene_data',
    'make_pubmed_gene_group',
    'write_boilerplate',
]

log = logging.getLogger(__name__)

# TODO make constants in :mod:`pybel.constants`
NAMESPACE_URL_FMT = 'DEFINE NAMESPACE {} AS URL "{}"'
NAMESPACE_PATTERN_FMT = 'DEFINE NAMESPACE {} AS PATTERN "{}"'
ANNOTATION_URL_FMT = 'DEFINE ANNOTATION {} AS URL "{}"'
ANNOTATION_PATTERN_FMT = 'DEFINE ANNOTATION {} AS PATTERN "{}"'

#: Allows for querying the Entrez Gene Summary utility by formatting with an entrez id or list of comma seperated ids
PUBMED_GENE_QUERY_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={}'


def split_document(lines):
    """Splits the lines over a document into the documents, definitions, and statements section

    :param iter[str] lines: A file or file-like that is an iterable over the lines of a document
    """
    lines = list(lines)
    end_document_section = 1 + max(i for i, line in enumerate(lines) if line.startswith('SET DOCUMENT'))
    end_definitions_section = 1 + max(i for i, line in enumerate(lines) if
                                      line.startswith('DEFINE ANNOTATION') or line.startswith('DEFINE NAMESPACE'))
    documents = [line for line in islice(lines, end_document_section) if not line.startswith('#')]
    definitions = [line for line in islice(lines, end_document_section, end_definitions_section) if
                   not line.startswith('#')]

    statements = lines[end_definitions_section:]

    return documents, definitions, statements


def merge(output_path, input_paths, merged_name=None, merged_contact=None, merged_description=None, merged_author=None):
    """Merges multiple BEL documents and maintains author information in comments

    Steps:

    1. Load all documents
    2. Identify document metadata information and ns/annot defs
    3. Postpend all statement groups with "- {author email}" and add comments with document information

    :param str output_path: Path to file to write merged BEL document
    :param iter[str] input_paths: List of paths to input BEL document files
    :param str merged_name: name for combined document
    :param str merged_contact: contact information for combine document
    :param str merged_description: description of combine document
    """
    metadata, defs, statements = [], [], []

    for input_path in input_paths:
        with open(os.path.expanduser(input_path)) as file:
            a, b, c = split_document(line.strip() for line in file)
            metadata.append(a)
            defs.append(set(b))
            statements.append(c)

    merged_contact = merged_contact if merged_contact is not None else ''
    merged_name = merged_name if merged_name is not None else 'MERGED DOCUMENT'
    merged_description = merged_description if merged_description is not None else 'This is a merged document'
    merged_author = merged_author if merged_author is not None else ''

    with open(os.path.expanduser(output_path), 'w') as file:
        for line in make_document_metadata(merged_name, merged_contact, merged_description,
                                           merged_author):
            print(line, file=file)

        for line in sorted(set().union(*defs)):
            print(line, file=file)

        for md, st in zip(metadata, statements):
            print(file=file)

            for line in md:
                print('# SUBDOCUMENT {}'.format(line), file=file)

            print(file=file)

            for line in st:
                print(line, file=file)


def make_document_metadata(name, contact, description, authors, version=None, copyright=None, licenses=None):
    """Builds a list of lines for the document metadata section of a BEL document

    :param str name: The unique name for this BEL document
    :param str contact: The email address of the maintainer
    :param str description: A description of the contents of this document
    :param str authors: The authors of this document
    :param str version: The version. Defaults to date in ``YYYYMMDD`` format.
    :param str copyright: Copyright information about this document
    :param str licenses: The license applied to this document
    :return: An iterator over the lines for the document metadata section
    :rtype: iter[str]
    """
    yield '# This document was created by PyBEL v{} on {}\n'.format(get_version(), time.asctime())

    yield '#' * 80
    yield '# Metadata'
    yield '#' * 80 + '\n'

    yield 'SET DOCUMENT Name = "{}"'.format(name)
    yield 'SET DOCUMENT Description = "{}"'.format(description.replace('\n', ''))
    yield 'SET DOCUMENT Version = "{}"'.format(version if version else time.strftime('%Y%m%d'))
    yield 'SET DOCUMENT Authors = "{}"'.format(authors)
    yield 'SET DOCUMENT ContactInfo = "{}"'.format(contact)

    if licenses:
        yield 'SET DOCUMENT Licenses = "{}"'.format(licenses)

    if copyright:
        yield 'SET DOCUMENT Copyright = "{}"'.format(copyright)

    yield ''


def make_document_namespaces(namespace_dict=None, namespace_patterns=None):
    """Builds a list of lines for the namespace definitions

    :param dict[str,str] namespace_dict: dictionary of {str name: str URL} of namespaces
    :param dict[str,str] namespace_patterns: A dictionary of {str name: str regex}
    :return: An iterator over the lines for the namespace definitions
    :rtype: iter[str]
    """
    namespace_dict = default_namespaces if namespace_dict is None else namespace_dict
    namespace_patterns = default_namespace_patterns if namespace_patterns is None else namespace_patterns

    yield '#' * 80
    yield '# Namespaces'
    yield '#' * 80 + '\n'
    yield '# Enumerated Namespaces\n'

    for name, url in sorted(namespace_dict.items()):
        yield NAMESPACE_URL_FMT.format(name, url)

    if namespace_patterns:
        yield '\n# Regular Expression Namespaces\n'

        for name, pattern in sorted(namespace_patterns.items()):
            yield NAMESPACE_PATTERN_FMT.format(name, pattern)

    yield ''


def make_document_annotations(annotation_dict=None, annotation_patterns=None):
    """Builds a list of lines for the annotation definitions

    :param dict[str,str] annotation_dict: A dictionary of {str name: str URL} of annotations
    :param dict[str,str] annotation_patterns: A dictionary of {str name: str regex}
    :return: An iterator over the lines for the annotation definitions
    :rtype: iter[str]
    """
    annotation_dict = default_annotations if annotation_dict is None else annotation_dict

    yield '#' * 80
    yield '# Annotations'
    yield '#' * 80 + '\n'

    for name, url in sorted(annotation_dict.items()):
        yield ANNOTATION_URL_FMT.format(name, url)

    if annotation_patterns:
        for name, pattern in sorted(annotation_patterns.items()):
            yield ANNOTATION_PATTERN_FMT.format(name, pattern)

    yield ''


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


def sanitize(s):
    if s is None:
        return None
    return s.strip().replace('\n', '')


def get_entrez_gene_data(entrez_ids):
    """Gets gene info from Entrez"""
    url = PUBMED_GENE_QUERY_URL.format(','.join(str(x).strip() for x in entrez_ids))
    response = requests.get(url)
    tree = ElementTree.fromstring(response.content)

    return {
        x.attrib['uid']: {
            'summary': sanitize(x.find('Summary').text),
            'description': x.find('Description').text
        }
        for x in tree.findall('./DocumentSummarySet/DocumentSummary')
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


def write_boilerplate(document_name, contact, description, authors, version=None, copyright=None,
                      licenses=None, namespace_dict=None, namespace_patterns=None, annotations_dict=None,
                      annotations_patterns=None, pmids=None, entrez_ids=None, file=None):
    """Writes a boilerplate BEL document, with standard document metadata, definitions. Optionally, if a
    list of PubMed identifiers are given, the citations and abstracts will be written for each.

    :param str document_name: The unique name for this BEL document
    :param str contact: The email address of the maintainer
    :param str description: A description of the contents of this document
    :param str authors: The authors of this document
    :param str version: The version. Defaults to current date in format YYYYMMDD.
    :param str copyright: Copyright information about this document
    :param str licenses: The license applied to this document
    
    :param dict[str, str] namespace_dict: an optional dictionary of {str name: str URL} of namespaces
    :param dict[str, str] namespace_patterns: An optional dictionary of {str name: str regex} namespaces
    :param dict[str, str] annotations_dict: An optional dictionary of {str name: str URL} of annotations
    :param dict[str, str] annotations_patterns: An optional dictionary of {str name: str regex} annotations
    :param iter[str] or iter[int] pmids: A list of PubMed identifiers to auto-populate with citation and abstract
    :param iter[str] or iter[int] entrez_ids: A list of Entrez identifiers to autopopulate the gene summary as evidence
    :param file file: A writable file or file-like. If None, defaults to :data:`sys.stdout`
    """
    file = sys.stdout if file is None else file

    metadata_iter = make_document_metadata(
        name=document_name,
        contact=contact,
        description=description,
        authors=authors,
        version=version,
        copyright=copyright,
        licenses=licenses
    )

    for line in metadata_iter:
        print(line, file=file)

    for line in make_document_namespaces(namespace_dict, namespace_patterns=namespace_patterns):
        print(line, file=file)

    for line in make_document_annotations(annotations_dict, annotation_patterns=annotations_patterns):
        print(line, file=file)

    print('#' * 80, file=file)
    print('# Statements', file=file)
    print('#' * 80 + '\n', file=file)

    if pmids is not None:
        for line in make_pubmed_abstract_group(pmids):
            print(line, file=file)

    if entrez_ids is not None:
        for line in make_pubmed_gene_group(entrez_ids):
            print(line, file=file)


def lint_helper(in_file, out_file=None):
    """Helps remove extraneous whitespace from the lines of a file

    :param file in_file: A readable file or file-like
    :param file out_file: A writable file or file-like
    """
    for line in in_file:
        print(line.strip(), file=out_file)


def lint_directory(source, target):
    """Adds a linted version of each document in the source directory to the target directory

    :param str source: Path to directory to lint
    :param str target: Path to directory to output
    """
    for path in os.listdir(source):
        if not path.endswith('.bel'):
            continue

        log.info('linting: %s', path)
        with open(os.path.join(source, path)) as i, open(os.path.join(target, path), 'w') as o:
            lint_helper(i, o)
