# -*- coding: utf-8 -*-

from __future__ import print_function

import logging

from pybel.constants import *
from pybel.manager.citation_utils import get_citations_by_pmids
from pybel.parser.parse_control import set_citation_stub, set_tag

__all__ = (
    'PubMedAnnotator',
    'rewrite_path_with_citations',
)

log = logging.getLogger(__name__)


class PubMedAnnotator:
    """Wraps functionality for rewriting BEL scripts"""

    def __init__(self, manager=None, group_size=200):
        self.cache = {}
        self.parser = set_tag + set_citation_stub
        self.group_size = group_size
        self.manager = manager

    def get_citations(self, pubmed_identifiers):
        """

        :param iter[str] pubmed_identifiers:
        """

        result = {}

        to_query = []

        for pmid in pubmed_identifiers:
            if pmid in self.cache:
                result[pmid] = self.cache[pmid]
            else:
                to_query.append(pmid)

        tr = get_citations_by_pmids(to_query, manager=self.manager, group_size=self.group_size)

        logging.warning('len tr: %s', len(tr))

        self.cache.update(**tr)

        result.update(**tr)

        return result

    def rewrite(self, lines):
        """Rewrites a BEL document

        Uses a clever stack so only the last 200 PubMed identifiers are ever in the memory and all the lines in between

        :param iter[str] lines: An iterable over the lines of a BEL file
        :return: An iterable of the new lines of a BEL file
        :rtype: iter[str]
        """
        pmids = set()
        stack = []

        for line in lines:
            try:
                r = self.parser.parseString(line)

                if r[0] != CITATION_TYPE_PUBMED:
                    stack.append(line)

                else:
                    pmid = int(r[1 if len(r) == 2 else 2])
                    log.warning('pmid: %s', pmid)
                    pmids.add(pmid)
                    stack.append(pmid)

                    if len(pmids) == self.group_size:
                        for item in self.help_iterate(stack, pmids):
                            yield item
            except:
                stack.append(line)

        for item in self.help_iterate(stack, pmids):
            yield item

    def help_iterate(self, stack, pmids):
        results = self.get_citations(pmids)

        for line in stack:
            if isinstance(line, int):
                pmid = str(line)

                result = results[pmid]

                yield 'SET Citation = {{"{}", "{}", "{}", "{}", "{}"}}'.format(
                    CITATION_TYPE_PUBMED,
                    result['title'],
                    pmid,
                    result[CITATION_DATE],
                    '|'.join(result[CITATION_AUTHORS]),
                )
            else:
                yield line

        pmids.clear()
        stack.clear()


def rewrite_path_with_citations(path, annotator=None, file=None):
    """Rewrites the BEL file at the given path with better citations

    :param str path:
    :param PubMedAnnotator annotator:
    :param file file:
    :return:
    """
    annotator = PubMedAnnotator() if annotator is None else annotator

    with open(path) as f:
        for line in annotator.rewrite(f):
            print(line.strip(), file=file)
