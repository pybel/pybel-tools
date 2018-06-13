# -*- coding: utf-8 -*-

"""Additional input and output methods"""

import logging
import os

from pybel import from_path, from_pickle, to_pickle, union
from pybel.manager import Manager

__all__ = [
    'from_path_ensure_pickle',
    'from_directory',
    'from_directory_pickles',
    'get_corresponding_gpickle_path',
]

log = logging.getLogger(__name__)

_bel_extension = '.bel'
_gpickle_extension = '.gpickle'
_json_extension = '.json'


def get_corresponding_gpickle_path(path):
    return path[:-len(_bel_extension)] + _gpickle_extension


def get_corresponding_json_path(path):
    return path[:-len(_bel_extension)] + _json_extension


def from_path_ensure_pickle(path, connection=None, **kwargs):
    """Parses a path exactly like :func:`pybel.from_path` unless a corresponding .gpickle file is available

    :param str path: A file path
    :param connection: database connection string to cache, pre-built :class:`Manager`, or None to use default cache
    :type connection: Optional[str or pybel.manager.Manager]
    :param kwargs:
    :rtype: pybel.BELGraph
    """
    if not path.endswith(_bel_extension):
        raise ValueError('Wrong extension. Should be .bel for file: {}'.format(path))

    gpickle_path = get_corresponding_gpickle_path(path)

    if os.path.exists(gpickle_path):
        with open(gpickle_path, 'rb') as file:
            return from_pickle(file)

    manager = Manager.ensure(connection=connection)
    graph = from_path(path, manager=manager, **kwargs)

    with open(gpickle_path, 'wb') as file:
        to_pickle(graph, file=file)

    return graph


def iter_pickle_paths_from_directory(directory, blacklist=None):
    """Iterates over file paths in a directory that are gpickles

    :param str directory: The directory
    :param Optional[list[str]] blacklist: An optional list of file names not to use
    :rtype: iter[str]
    :raises: ValueError
    """
    if not os.path.isdir(directory):
        raise ValueError('Not a directory: {}'.format(directory))

    for path in os.listdir(directory):
        if blacklist and path in blacklist:
            continue

        if not path.endswith('.gpickle'):
            continue

        yield os.path.join(directory, path)


def iter_from_pickles(paths):
    """Iterates over the pickled BEL graphs in a directory

    :param iter[str] paths:
    :rtype: iter[pybel.BELGraph]
    """
    for path in paths:
        if not path.endswith('.gpickle'):
            log.info('Wrong extension. Should be .gpickle: %s', path)
            continue
        yield from_pickle(path)


def iter_from_pickles_from_directory(directory, blacklist=None):
    """Iterates over the pickled BEL graphs in a directory

    :param str directory:
    :param Optional[list[str]] blacklist: An optional list of file names not to use
    :rtype: iter[pybel.BELGraph]
    """
    return iter_from_pickles(iter_pickle_paths_from_directory(directory, blacklist=blacklist))


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


def iter_paths_from_directory(directory):
    for filename in os.listdir(directory):
        if not filename.endswith(_bel_extension):
            continue
        yield filename


def iter_from_directory(directory, connection=None):
    """Parses all BEL scripts in the given directory with :func:`load_paths` and returns the union of the resulting
    graphs.

    :param str directory: A path to a directory
    :param connection: database connection string to cache, pre-built :class:`Manager`, or None to use default cache
    :type connection: Optional[str or pybel.manager.Manager]
    :rtype: iter[pybel.BELGraph]
    """
    for filename in iter_paths_from_directory(directory):
        path = os.path.join(directory, filename)
        yield from_path_ensure_pickle(path, connection=connection)


def from_directory(directory, connection=None):
    """Parses all BEL scripts in the given directory with :func:`load_paths` and returns the union of the resulting
    graphs.

    :param str directory: A path to a directory
    :param connection: database connection string to cache, pre-built :class:`Manager`, or None to use default cache
    :type connection: Optional[str or pybel.manager.Manager]
    :rtype: pybel.BELGraph
    """
    return union(iter_from_directory(directory, connection=connection))
