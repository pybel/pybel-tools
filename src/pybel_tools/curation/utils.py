# -*- coding: utf-8 -*-

"""Curation utilities."""

from typing import Iterable, List, Optional

import click
from gitlab.v4.objects import Project

from pybel.manager.citation_utils import get_pubmed_citation_response

min_year_option = click.option(
    '--min-year', type=int, default=2001, show_default=True,
    help='Minimum publication year. Before 2001 is dubious.',
)


def make_issues_from_pmids(
    project: Project,
    pmids: Iterable[str],
    min_year: Optional[int] = None,
    allow_closed: bool = False,
    make_issues: bool = False,
    labels: Optional[List[str]] = None,
) -> None:
    """Make issues on the GitLab project for the given articles."""
    results = get_pubmed_citation_response(pmids)['result']
    existing_issue_titles = {
        issue.title
        for issue in project.issues.list(all=True)
    }

    for pmid in pmids:
        result = results[pmid]
        first_author_surname = result['sortfirstauthor'].split()[0].lower()
        pubyear = int(result['sortpubdate'].split('/')[0])
        if min_year is not None and pubyear < min_year:
            continue

        pmc = get_pmc_from_result(result)
        if pmc is not None:
            issue_title = f'''{first_author_surname}{pubyear} (pmid:{pmid}, pmc:{pmc}) "{result['title']}"'''
        elif allow_closed:
            issue_title = f'''{first_author_surname}{pubyear} (pmid:{pmid}) "{result['title']}"'''
        else:
            continue

        if make_issues and issue_title not in existing_issue_titles:
            project_issue = project.issues.create({
                'title': f'Curate {issue_title}',
                # 'description': 'Something useful here.'
            })
            if pmc is not None:
                project_issue.labels.append('Availability: PMC')
            project_issue.labels = ['Curation']
            if labels is not None:
                project_issue.labels.extend(labels)
            project_issue.save()


def get_pmc_from_result(result) -> Optional[str]:
    for article_id in result['articleids']:
        if article_id['idtype'] == 'pmc':
            return article_id['value']


PDF_PREFIX = 'PDF:'
DOWNLOAD_PREFIX = 'Download:'


def get_issue_pdf(description: str) -> str:
    """Get the link to the PDF for a given issue."""
    for line in description.splitlines():
        line = lstrip_list(line)
        if line.startswith(PDF_PREFIX) and len(line) > len(PDF_PREFIX):
            return line[len(PDF_PREFIX):]
        if line.startswith(DOWNLOAD_PREFIX) and len(line) > len(DOWNLOAD_PREFIX):
            return line[len(DOWNLOAD_PREFIX):]


def lstrip_list(s: str) -> str:
    """Left strip a string of its bullet point."""
    return s.strip().lstrip('-').lstrip()
