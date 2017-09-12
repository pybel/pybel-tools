# -*- coding: utf-8 -*-

import logging
import time
from datetime import datetime

import re
import requests

from pybel.constants import *
from pybel.manager import Manager
from pybel.parser.parse_control import set_citation_stub, set_tag
from pybel_tools.utils import grouper

__all__ = (
    'get_citations_by_pmids',
)

log = logging.getLogger(__name__)

EUTILS_URL_FMT = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id={}"

re1 = re.compile('^[12][0-9]{3} [a-zA-Z]{3} \d{1,2}$')
re2 = re.compile('^[12][0-9]{3} [a-zA-Z]{3}$')
re3 = re.compile('^[12][0-9]{3}$')
re4 = re.compile('^[12][0-9]{3} [a-zA-Z]{3}-[a-zA-Z]{3}$')


def get_citations_by_pmids(pmids, group_size=200, sleep_time=1, return_errors=False, manager=None):
    """Gets the citation information for the given list of PubMed identifiers using the NCBI's eutils service

    :param iter[str] or iter[int] pmids: an iterable of PubMed identifiers
    :param int group_size: The number of PubMed identifiers to query at a time
    :param int sleep_time: Number of seconds to sleep between queries
    :param bool return_errors: Should a set of erroneous PubMed identifiers be returned?
    :param manager: An RFC-1738 database connection string, a pre-built :class:`pybel.manager.Manager`,
                    or ``None`` for default connection
    :type manager: None or str or Manager
    :return: A dictionary of {pmid: pmid data dictionary} or a pair of this dictionary and a set ot erroneous
            pmids if return_errors is :data:`True`
    :rtype: dict[str,dict]
    """
    pmids = {str(pmid).strip() for pmid in pmids}
    log.info('querying %d PubMed identifiers', len(pmids))

    manager = Manager.ensure(manager)

    result = {}
    unresolved_pmids = {}

    for pmid in pmids:
        citation = manager.get_or_create_citation(type=CITATION_TYPE_PUBMED, reference=pmid)

        if not citation.date or not citation.name or not citation.authors:
            unresolved_pmids[pmid] = citation
            continue

        result[pmid] = citation.to_json()

    manager.session.commit()

    log.info('Used %d citations from cache', len(pmids) - len(unresolved_pmids))

    if not unresolved_pmids:
        return (result, set()) if return_errors else result

    log.info('Looking up %d citations in PubMed', len(unresolved_pmids))

    errors = set()
    t = time.time()

    for pmid_list in grouper(group_size, unresolved_pmids):
        url = EUTILS_URL_FMT.format(','.join(pmid for pmid in pmid_list if pmid))
        response = requests.get(url).json()

        for pmid in response['result']['uids']:
            p = response['result'][pmid]

            if 'error' in p:
                log.warning("Error downloading PubMed identifier: %s", pmid)
                errors.add(pmid)
                continue

            result[pmid] = {
                'title': p['title'],
                'last': p['lastauthor'],
                CITATION_NAME: p['fulljournalname'],
                'volume': p['volume'],
                'issue': p['issue'],
                'pages': p['pages'],
                'first': p['sortfirstauthor'],
            }

            citation = unresolved_pmids[pmid]

            citation.name = result[pmid][CITATION_NAME]
            citation.title = result[pmid]['title']
            citation.volume = result[pmid]['volume']
            citation.issue = result[pmid]['issue']
            citation.pages = result[pmid]['pages']
            citation.first = result[pmid]['first']
            citation.last = result[pmid]['last']

            if 'authors' in p:
                result[pmid][CITATION_AUTHORS] = [author['name'] for author in p['authors']]
                for author in result[pmid][CITATION_AUTHORS]:
                    author_model = manager.get_or_create_author(author)
                    if author_model not in citation.authors:
                        citation.authors.append(author_model)

            publication_date = p['pubdate']

            if re1.search(publication_date):
                result[pmid][CITATION_DATE] = datetime.strptime(p['pubdate'], '%Y %b %d').strftime('%Y-%m-%d')
            elif re2.search(publication_date):
                result[pmid][CITATION_DATE] = datetime.strptime(p['pubdate'], '%Y %b').strftime('%Y-%m-01')
            elif re3.search(publication_date):
                result[pmid][CITATION_DATE] = p['pubdate'] + "-01-01"
            elif re4.search(publication_date):
                result[pmid][CITATION_DATE] = datetime.strptime(p['pubdate'][:-4], '%Y %b').strftime('%Y-%m-01')
            else:
                log.info('Date with strange format: %s', p['pubdate'])

            if CITATION_DATE in result[pmid]:
                citation.date = datetime.strptime(result[pmid][CITATION_DATE], '%Y-%m-%d')

        manager.session.commit()  # commit in groups

        # Don't want to hit that rate limit
        time.sleep(sleep_time)

    log.info('retrieved %d PubMed identifiers in %.02f seconds', len(unresolved_pmids), time.time() - t)

    return (result, errors) if return_errors else result


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
                        yield from self.help_iterate(stack, pmids)
            except:
                stack.append(line)

        yield from self.help_iterate(stack, pmids)

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
