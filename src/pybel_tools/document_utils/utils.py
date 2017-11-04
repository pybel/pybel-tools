# -*- coding: utf-8 -*-

import logging
import os

__all__ = [
    'lint_file',
    'lint_directory',
]

log = logging.getLogger(__name__)


def lint_file(in_file, out_file=None):
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
            lint_file(i, o)
