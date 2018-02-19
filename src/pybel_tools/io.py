# -*- coding: utf-8 -*-


import logging
import os

from pybel import Manager, from_path, from_pickle, union

__all__ = [
    'from_directory',
    'from_directory_pickles',
]

log = logging.getLogger(__name__)


def iter_pickle_paths_from_directory(directory):
    """Iterates over file paths in a directory that are gpickles

    :param str directory: The directory
    :rtype: iter[str]
    :raises: ValueError
    """
    if not os.path.isdir(directory):
        raise ValueError('Not a directory: {}'.format(directory))

    for path in os.listdir(directory):
        if path.endswith('.gpickle'):
            yield path


def iter_from_pickles(paths):
    """Iterates over the pickled BEL graphs in a directory

    :param iter[str] paths:
    :rtype: iter[pybel.BELGraph]
    """
    for path in paths:
        if not path.endswith('.gpickle'):
            log.info('not a gpickle: %s', path)
            continue
        yield from_pickle(path)


def iter_from_pickles_from_directory(directory):
    """Iterates over the pickled BEL graphs in a directory

    :param str directory:
    :rtype: iter[pybel.BELGraph]
    """
    return iter_from_pickles(iter_pickle_paths_from_directory(directory))


def from_pickles(paths):
    """Loads multiple PyBEL pickles with :func:`pybel.from_pickle` and returns the union of the resulting graphs.

    :param iter[str] paths: An iterable over paths to PyBEL pickles
    :rtype: pybel.BELGraph
    """
    return union(iter_from_pickles(paths))


def from_directory_pickles(directory):
    """Parses all BEL scripts in the given directory with :func:`load_paths` and returns the union of the resulting
    graphs.

    :param str directory: A path to a directory
    :rtype: pybel.BELGraph
    """
    return union(iter_from_pickles_from_directory(directory))


def from_directory(directory, connection=None):
    """Parses all BEL scripts in the given directory with :func:`load_paths` and returns the union of the resulting
    graphs.

    :param str directory: A path to a directory
    :param connection: A custom database connection string or manager
    :type connection: Optional[str or pybel.manager.Manager]
    :rtype: pybel.BELGraph
    """
    paths = [
        path
        for path in os.listdir(directory)
        if path.endswith('.bel')
    ]
    log.info('loadings %d paths: %s', len(paths), ', '.join(paths))
    return load_paths(paths, connection=connection)


def load_paths(paths, connection=None):
    """Parses multiple BEL scripts with :func:`pybel.from_path` and returns the union of the resulting graphs.

    :param iter[str] paths: An iterable over paths to BEL scripts
    :param connection: A custom database connection string or manager
    :type connection: Optional[str or pybel.manager.Manager]
    :rtype: pybel.BELGraph
    """
    manager = Manager.ensure(connection)

    return union(
        from_path(path, manager=manager)
        for path in paths
    )
