# -*- coding: utf-8 -*-

"""Utilities for loading and exporting BEL graphs"""

import logging

import os
import requests

import pybel
from pybel import from_path, to_pickle, from_pickle, to_database
from pybel.io.io_exceptions import ImportVersionWarning
from pybel.manager import Manager
from pybel.struct import union
from pybel.utils import get_version as get_pybel_version
from .constants import DEFAULT_SERVICE_URL
from .integration import HGNCAnnotator
from .mutation import opening_on_central_dogma
from .mutation.metadata import enrich_pubmed_citations
from .selection import get_subgraph_by_annotation_value
from .summary import get_annotation_values

__all__ = [
    'load_paths',
    'load_directory',
    'subgraphs_to_pickles',
    'convert_paths',
    'convert_directory',
    'to_pybel_web',
]

log = logging.getLogger(__name__)


def load_paths(paths, connection=None):
    """Loads a group of BEL graphs.

    Internally, this function uses a shared :class:`pybel.parser.MetadataParser` to cache the definitions more
    efficiently.

    :param iter[str] paths: An iterable over paths to BEL scripts
    :param str connection: A custom database connection string
    :type connection: None or str or pybel.manager.Manager
    :return: A BEL graph comprised of the union of all BEL graphs produced by each BEL script
    :rtype: pybel.BELGraph
    """
    manager = Manager.ensure(connection)

    return union(
        from_path(path, manager=manager)
        for path in paths
    )


def load_directory(directory, connection=None):
    """Compiles all BEL scripts in the given directory and returns as a merged BEL graph using :func:`load_paths`

    :param str directory: A path to a directory
    :param str connection: A custom database connection string
    :type connection: None or str or pybel.manager.Manager
    :return: A BEL graph comprised of the union of all BEL graphs produced by each BEL script
    :rtype: pybel.BELGraph
    """
    paths = (
        path
        for path in os.listdir(directory)
        if path.endswith('.bel')
    )

    return load_paths(paths, connection=connection)


def get_paths_recursive(directory, extension='.bel', exclude_directory_pattern=None):
    """Gets all file paths in a given directory to BEL documents

    :param str directory: The base directory to walk
    :param str extension: Extensions of files to keep
    :param str exclude_directory_pattern: Any directory names to exclude
    """
    for root, directory, files in os.walk(directory):
        if exclude_directory_pattern and exclude_directory_pattern in directory:
            continue

        for file in files:
            if file.endswith(extension):
                yield os.path.join(root, file)


def convert_paths(paths, connection=None, upload=False, pickle=False, store_parts=False,
                  infer_central_dogma=True, enrich_citations=False, enrich_genes=False, enrich_go=False, send=False,
                  version_in_path=False, **kwargs):
    """Recursively parses and either uploads/pickles graphs in a given set of files

    :param iter[str] paths: The paths to convert
    :param connection: The connection
    :type connection: None or str or pybel.manager.Manager
    :param bool upload: Should the networks be uploaded to the cache?
    :param bool pickle: Should the networks be saved as pickles?
    :param bool store_parts: Should the networks be uploaded to the edge store?
    :param bool infer_central_dogma: Should the central dogma be inferred for all proteins, RNAs, and miRNAs
    :param bool enrich_citations: Should the citations be enriched using Entrez Utils?
    :param bool enrich_genes: Should the genes' descriptions be downloaded from Gene Cards?
    :param bool enrich_go: Should the biological processes' descriptions be downloaded from Gene Ontology?
    :param bool send: Send to PyBEL Web?
    :param bool version_in_path: Add the current pybel version to the pathname
    :param kwargs: Parameters to pass to :func:`pybel.from_path`
    """
    from .integration.description.go_annotator import GOAnnotator
    manager = Manager.ensure(connection)
    hgnc_annotator = HGNCAnnotator(preload=enrich_genes)
    go_annotator = GOAnnotator(preload=enrich_go)

    failures = []

    for path in paths:
        log.info('Parsing path: %s', path)

        try:
            network = from_path(path, manager=manager, **kwargs)
        except Exception as e:
            log.exception('Problem parsing %s', path)
            failures.append((path, e))
            continue

        if infer_central_dogma:
            opening_on_central_dogma(network)

        if enrich_citations:
            enrich_pubmed_citations(network, manager=manager)

        if enrich_genes and hgnc_annotator.download_successful:
            hgnc_annotator.annotate(network)

        if enrich_go and go_annotator.download_successful:
            go_annotator.annotate(network)

        if upload or store_parts:
            to_database(network, connection=manager, store_parts=store_parts)

        if pickle:
            name = path[:-4]  # [:-4] gets rid of .bel at the end of the file name

            if version_in_path:
                new_path = '{}-{}.gpickle'.format(name, get_pybel_version())
            else:
                new_path = '{}.gpickle'.format(name)

            to_pickle(network, new_path)

        if send:
            to_pybel_web(network)

    return failures


def convert_directory(directory, connection=None, upload=False, pickle=False, store_parts=False,
                      infer_central_dogma=True, enrich_citations=False, enrich_genes=False, enrich_go=False, send=False,
                      exclude_directory_pattern=None, version_in_path=False, **kwargs):
    """Recursively parses and either uploads/pickles graphs in a given directory and sub-directories

    :param str directory: The directory to look through
    :param connection: The connection
    :type connection: None or str or pybel.manager.Manage.
    :param bool upload: Should the networks be uploaded to the cache?
    :param bool pickle: Should the networks be saved as pickles?
    :param bool store_parts: Should the networks be uploaded to the edge store?
    :param bool infer_central_dogma: Should the central dogma be inferred for all proteins, RNAs, and miRNAs
    :param bool enrich_citations: Should the citations be enriched using Entrez Utils?
    :param bool enrich_genes: Should the genes' descriptions be downloaded from Gene Cards?
    :param bool enrich_go: Should the biological processes' descriptions be downloaded from Gene Ontology?
    :param bool send: Send to PyBEL Web?
    :param str exclude_directory_pattern: A pattern to use to skip directories
    :param bool version_in_path: Add the current pybel version to the pathname
    :param kwargs: Parameters to pass to :func:`pybel.from_path`
    """
    paths = list(get_paths_recursive(directory, exclude_directory_pattern=exclude_directory_pattern))
    log.info('Paths to parse: %s', paths)

    result = convert_paths(
        paths,
        connection=connection,
        upload=upload,
        pickle=pickle,
        store_parts=store_parts,
        infer_central_dogma=infer_central_dogma,
        enrich_citations=enrich_citations,
        enrich_genes=enrich_genes,
        enrich_go=enrich_go,
        send=send,
        version_in_path=version_in_path,
        **kwargs
    )

    return result


def upload_recursive(directory, connection=None, store_parts=False, exclude_directory_pattern=None):
    """Recursively uploads all gpickles in a given directory and sub-directories
    
    :param str directory: the directory to traverse
    :param connection: A connection string or manager
    :type connection: None or str or pybel.manage.Manager
    :param bool store_parts: Should the edge store be used?
    :param str exclude_directory_pattern: Any directory names to exclude
    """
    manager = Manager.ensure(connection)
    paths = list(get_paths_recursive(
        directory,
        extension='.gpickle',
        exclude_directory_pattern=exclude_directory_pattern
    ))
    log.info('Paths to upload: %s', paths)

    for path in paths:
        try:
            network = from_pickle(path)
        except (ImportError, ImportVersionWarning):
            log.warning('%s uses a pickle from an old version of PyBEL. Quitting.', path)
            continue

        to_database(network, connection=manager, store_parts=store_parts)


def subgraphs_to_pickles(network, directory=None, annotation='Subgraph'):
    """Groups the given graph into subgraphs by the given annotation with :func:`get_subgraph_by_annotation` and
    outputs them as gpickle files to the given directory with :func:`pybel.to_pickle`

    :param pybel.BELGraph network: A BEL network
    :param str directory: A directory to output the pickles
    :param str annotation: An annotation to split by. Suggestion: ``Subgraph``
    """
    directory = os.getcwd() if directory is None else directory
    for value in get_annotation_values(network, annotation=annotation):
        sg = get_subgraph_by_annotation_value(network, annotation, value)
        sg.document.update(network.document)

        file_name = '{}_{}.gpickle'.format(annotation, value.replace(' ', '_'))
        path = os.path.join(directory, file_name)
        to_pickle(sg, path)


def to_pybel_web(network, service=None):
    """Sends a graph to the receiver service and returns the :mod:`requests` response object

    :param pybel.BELGraph network: A BEL network
    :param str service: The location of the PyBEL web server. Defaults to :data:`DEFAULT_SERVICE_URL`
    :return: The response object from :mod:`requests`
    :rtype: requests.Response
    """
    service = DEFAULT_SERVICE_URL if service is None else service
    url = service + '/api/receive'
    headers = {'content-type': 'application/json'}
    return requests.post(url, json=pybel.to_json(network), headers=headers)
