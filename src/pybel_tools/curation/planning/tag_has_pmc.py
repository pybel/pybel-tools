# -*- coding: utf-8 -*-

"""A one-off script for updating availability tags in GitLab.

This script gets all curation issues that are tagged as ``Availability: Missing``, looks up their
information in PubMed, and switches for ``Availability: PMC`` when appropriate.
"""

from typing import Any, Mapping, Optional

import click
from easy_config.contrib.click import args_from_config
from pybel.manager.citation_utils import get_pubmed_citation_response
from pybel_git.gitlab import GitlabConfig
from tqdm import tqdm

from ..utils import get_issue_pdf, get_pmc_from_result, lstrip_list

PMID_PREFIX = 'PMID:'
AVAILABILITY_PREFIX = 'Availability: '
AVAILABILE_PDF = 'Availability: PDF'
AVAILABILE_PMC = 'Availability: PMC'
AVAILABILE_MISSING = 'Availability: Missing'


@click.command()
@args_from_config(GitlabConfig)
def main(project_id: int, url: str, token: str):
    """Update the availability tags."""
    gitlab_config = GitlabConfig.load(  # noqa: S106
        project_id=project_id,
        url=url,
        token=token,
    )
    project = gitlab_config.get_project()
    issues = [
        issue
        for issue in project.issues.list(all=True)
        if 'Curation' in issue.labels
    ]

    missing_dict = {}
    for issue in issues:
        if get_issue_pdf(issue.description):
            if AVAILABILE_PDF not in issue.labels:
                issue.labels.append(AVAILABILE_PDF)
                issue.save()
            continue

        pmid = get_pmid_from_title(issue.title)

        if pmid is None:
            pmid = get_pmid_from_description(issue.description)

        if pmid is None:
            print(f'MISSING PMID: {issue.title}\n{issue.description}\n{"=" * 80}\n\n')
            continue

        availability = get_availability(issue)

        if availability is None:
            issue.labels.append(AVAILABILE_MISSING)
            issue.save()
            availability = 'Missing'

        if availability == 'Missing':
            if get_issue_pdf(issue.description):
                print(f'PDF available for (pmid:{pmid}) / {issue.title}')
                issue.labels.remove(AVAILABILE_MISSING)
                issue.labels.append(AVAILABILE_PDF)
                issue.save()
            else:
                print(f'Availability missing for (pmid:{pmid}) / {issue.title}')
                missing_dict[pmid] = issue

    print(f'{len(missing_dict)} PMIDs were missing availabilities')

    pmid_to_pmc = get_availabilities_from_eutils(missing_dict)

    print(f'Found {len(pmid_to_pmc)} in PMC.')
    it = tqdm(missing_dict.items(), desc='Updating issues')
    for pmid, issue in it:
        pmc = pmid_to_pmc.get(pmid)
        if pmc is None:
            continue
        it.write(f'Updating availability on {issue.title}')
        issue.labels.remove(AVAILABILE_MISSING)
        issue.labels.append(AVAILABILE_PMC)
        issue.discussions.create({
            'body': f'Available from PubMed Central at https://identifiers.org/pmc:{pmc}'
        })
        issue.save()


def get_availability(issue):
    for label in issue.labels:
        if label.startswith(AVAILABILITY_PREFIX):
            return label[len(AVAILABILITY_PREFIX):]


def get_pmid_from_title(title) -> Optional[str]:
    try:
        index = title.index('pmid:')
    except ValueError:
        return None
    else:
        if 0 <= index:
            return title[index + len('pmid:'):title.index(')', index)]


def get_pmid_from_description(description):
    for line in description.splitlines():
        line = lstrip_list(line).replace('[', '')
        if line.startswith(PMID_PREFIX):
            if ']' in line:
                return line[len(PMID_PREFIX) + 1:line.index(']')]
            else:
                return line[len(PMID_PREFIX):].strip()


def get_availabilities_from_eutils(pmids: Mapping[str, Any]) -> Mapping[str, str]:
    rv = {}
    results = get_pubmed_citation_response(pmids)['result']
    for pmid, issue in pmids.items():
        try:
            result = results[pmid]
        except KeyError:
            print(f'Error on {issue.title}')
            continue
        pmc = get_pmc_from_result(result)
        if pmc is not None:
            rv[pmid] = pmc
    return rv


if __name__ == '__main__':
    main()
