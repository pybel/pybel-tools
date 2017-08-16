# -*- coding: utf-8 -*-

import logging
import re
import time
from collections import defaultdict
from datetime import datetime

import requests

from pybel.constants import *
from pybel.manager import build_manager
from pybel.manager.models import Citation
from pybel.parser.parse_control import set_citation, set_tag

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
    :return: A dictionary of {pmid: pmid data dictionary} or a pair of this dictionary and a set ot erroneous
            pmids if return_errors is :data:`True`
    :rtype: dict[str,dict]
    """
    pmids = [str(pmid).strip() for pmid in sorted(pmids)]
    log.info('querying %d PubMed identifiers', len(pmids))
    result = defaultdict(dict)

    manager = build_manager(manager)

    unresolved_pmids = []

    for pmid in pmids:
        c = manager.session.query(Citation).filter(Citation.type == 'PubMed', Citation.reference == pmid).one_or_none()

        if c is None:
            unresolved_pmids.append(pmid)
            continue

        result[pmid] = {
            CITATION_AUTHORS: '|'.join(c.authors),
            CITATION_TYPE: CITATION_TYPE_PUBMED,
            CITATION_REFERENCE: pmid,
            CITATION_DATE: c.date,
            CITATION_NAME: c.name
        }

    errors = set()
    t = time.time()

    for pmidList in [','.join(unresolved_pmids[i:i + group_size]) for i in range(0, len(unresolved_pmids), group_size)]:
        url = EUTILS_URL_FMT.format(pmidList)
        response = requests.get(url).json()

        pmid_result = response['result']

        for pmid in pmid_result['uids']:
            p = pmid_result[pmid]

            if 'error' in p:
                log.warning("Error downloading PubMed identifier: %s", pmid)
                errors.add(pmid)
                continue

            if 'authors' in p:
                result[pmid][CITATION_AUTHORS] = [author['name'] for author in p['authors']]

            if re1.search(p['pubdate']):
                result[pmid][CITATION_DATE] = datetime.strptime(p['pubdate'], '%Y %b %d').strftime('%Y-%m-%d')
            elif re2.search(p['pubdate']):
                result[pmid][CITATION_DATE] = datetime.strptime(p['pubdate'], '%Y %b').strftime('%Y-%m-01')
            elif re3.search(p['pubdate']):
                result[pmid][CITATION_DATE] = p['pubdate'] + "-01-01"
            elif re4.search(p['pubdate']):
                result[pmid][CITATION_DATE] = datetime.strptime(p['pubdate'][:-4], '%Y %b').strftime('%Y-%m-01')
            else:
                log.info('Date with strange format: %s', p['pubdate'])

            result[pmid].update({
                'title': p['title'],
                'lastauthor': p['lastauthor'],
                CITATION_NAME: p['fulljournalname'],
                'volume': p['volume'],
                'issue': p['issue'],
                'pages': p['pages'],
                'firstauthor': p['sortfirstauthor'],
            })

            c = Citation(
                authors=[
                    manager.get_or_create_author(author)
                    for author in result[pmid][CITATION_AUTHORS]
                ],
                **{
                    CITATION_TYPE: CITATION_TYPE_PUBMED,
                    CITATION_NAME: result[pmid][CITATION_NAME],
                    CITATION_REFERENCE: pmid,
                    CITATION_DATE: datetime.strptime(result[pmid][CITATION_DATE], '%Y-%m-%d'),
                }
            )

            manager.session.add(c)
            manager.session.commit()

        # Don't want to hit that rate limit
        time.sleep(sleep_time)

    log.info('retrieved PubMed identifiers in %.02f seconds', time.time() - t)

    return (result, errors) if return_errors else result


class PubMedAnnotator:
    def __init__(self):
        self.cache = {}
        self.parser = set_tag + set_citation

    def get_citations(self, pubmed_identifiers):
        """

        :param iter[str] pubmed_identifiers:
        :rtype: iter[
        """

        result = {}

        to_query = []

        for pmid in pubmed_identifiers:
            if pmid in self.cache:
                result[pmid] = self.cache[pmid]
            else:
                to_query.append(pmid)

        tr = get_citations_by_pmids(to_query)

        logging.warning('len tr: %s', len(tr))

        self.cache.update(**tr)

        result.update(**tr)

        return result

    def rewrite(self, lines, group_size=200):
        """Rewrites a BEL document

        Uses a clever stack so only the last 200 pubmeds are ever in the memory and all the lines in between

        :param iter[str] lines: An iterable over the lines of a BEL file
        :param int group_size: The number of PubMed identifiers to query at a time
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

                    if len(pmids) == group_size:
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
