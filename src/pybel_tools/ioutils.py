# -*- coding: utf-8 -*-

"""Utilities for loading and exporting BEL graphs"""

import logging
import os

from pybel import from_pickle, to_database, to_pickle, to_web
from pybel.io.exc import ImportVersionWarning
from pybel.manager import Manager
from .io import from_path_ensure_pickle
from .mutation import add_canonical_names, enrich_pubmed_citations, infer_central_dogma as infer_central_dogma_mutator
from .selection import get_subgraph_by_annotation_value
from .summary import get_annotation_values

__all__ = [
    'subgraphs_to_pickles',
    'convert_paths',
    'convert_directory',
    'get_paths_recursive',
    'upload_recursive',
]

log = logging.getLogger(__name__)


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


def convert_paths(paths, connection=None, upload=False, canonicalize=True, infer_central_dogma=True,
                  enrich_citations=False, send=False, **kwargs):
    """Recursively parses and either uploads/pickles graphs in a given set of files

    :param iter[str] paths: The paths to convert
    :param connection: The connection
    :type connection: None or str or pybel.manager.Manager
    :param bool upload: Should the networks be uploaded to the cache?
    :param bool canonicalize: Calculate canonical nodes?
    :param bool infer_central_dogma: Should the central dogma be inferred for all proteins, RNAs, and miRNAs
    :param bool enrich_citations: Should the citations be enriched using Entrez Utils?
    :param bool send: Send to PyBEL Web?
    :param kwargs: Parameters to pass to :func:`pybel.from_path`
    :return: A pair of a dictionary {path: bel graph} and list of failed paths
    :rtype: tuple[dict[str,pybel.BELGraph],list[str]]
    """
    manager = Manager.ensure(connection)

    successes = {}
    failures = []

    for path in paths:
        try:
            graph = from_path_ensure_pickle(path, connection=manager, **kwargs)
        except Exception as e:
            log.exception('problem parsing %s', path)
            failures.append((path, e))
            continue

        if canonicalize:
            add_canonical_names(graph)

        if infer_central_dogma:
            infer_central_dogma_mutator(graph)

        if enrich_citations:
            enrich_pubmed_citations(graph=graph, manager=manager)

        if upload:
            to_database(graph, connection=manager, store_parts=True)

        if send:
            response = to_web(graph)
            log.info('sent to PyBEL Web with response: %s', response.json())

        successes[path] = graph

    return successes, failures


def convert_directory(directory, connection=None, upload=False, pickle=False, canonicalize=True,
                      infer_central_dogma=True, enrich_citations=False, enrich_genes=False, enrich_go=False, send=False,
                      exclude_directory_pattern=None, version_in_path=False, **kwargs):
    """Recursively parses and either uploads/pickles graphs in a given directory and sub-directories

    :param str directory: The directory to look through
    :param connection: The connection
    :type connection: None or str or pybel.manager.Manage.
    :param bool upload: Should the networks be uploaded to the cache?
    :param bool pickle: Should the networks be saved as pickles?
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
        canonicalize=canonicalize,
        infer_central_dogma=infer_central_dogma,
        enrich_citations=enrich_citations,
        enrich_genes=enrich_genes,
        enrich_go=enrich_go,
        send=send,
        version_in_path=version_in_path,
        **kwargs
    )

    return result


def upload_recursive(directory, connection=None, exclude_directory_pattern=None):
    """Recursively uploads all gpickles in a given directory and sub-directories
    
    :param str directory: the directory to traverse
    :param connection: A connection string or manager
    :type connection: Optional[str or pybel.manage.Manager]
    :param Optional[str] exclude_directory_pattern: Any directory names to exclude
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

        to_database(network, connection=manager, store_parts=True)


def subgraphs_to_pickles(network, annotation, directory=None):
    """Groups the given graph into subgraphs by the given annotation with :func:`get_subgraph_by_annotation` and
    outputs them as gpickle files to the given directory with :func:`pybel.to_pickle`

    :param pybel.BELGraph network: A BEL network
    :param str annotation: An annotation to split by. Suggestion: ``Subgraph``
    :param Optional[str] directory: A directory to output the pickles
    """
    directory = directory or os.getcwd()

    for value in get_annotation_values(network, annotation):
        sg = get_subgraph_by_annotation_value(network, annotation, value)
        sg.document.update(network.document)

        file_name = '{}_{}.gpickle'.format(annotation, value.replace(' ', '_'))
        path = os.path.join(directory, file_name)
        to_pickle(sg, path)
