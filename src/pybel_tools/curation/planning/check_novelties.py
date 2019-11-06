# -*- coding: utf-8 -*-

"""This script assesses the novelty of pending curation tasks.

Currently, is limited to articles where PMC is available to ensure
good INDRA coverage.
"""

import json
import logging
from typing import Optional

import click
from easy_config.contrib.click import args_from_config
from gitlab.v4.objects import Issue, Project
from hbp_knowledge import get_graph
from pybel_git.gitlab import GitlabConfig

from pybel import BELGraph
from pybel_tools.assess_completeness import CompletenessSummary, assess_completeness
from ..recuration.utils import CURATION_LABEL

_prefix = '- PMID: ['


@click.command()
@args_from_config(GitlabConfig)
@click.option('-o', '--output', type=click.File('w'))
def main(project_id: int, url: str, token: str, output) -> None:
    """Assess the completeness of HBP curation tasks with respect to CONIB."""
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('hbp').setLevel(logging.INFO)

    gitlab_config = GitlabConfig.load(  # noqa: S106
        project_id=project_id,
        url=url,
        token=token,
    )
    project = gitlab_config.get_project()
    do_it(project, output)


def do_it(project: Project, output):
    graph = get_graph()

    summaries = assess_project_completeness(project=project, graph=graph)

    if output is not None:
        json.dump(list(summaries), output, indent=2)
    else:
        for summary in summaries:
            click.echo(json.dumps(summary, indent=2))


def assess_project_completeness(*, project: Project, graph: BELGraph):
    """Summarize thee novelty of all issues in the project."""
    issues = project.issues.list(labels=[CURATION_LABEL])
    for issue in issues:
        click.echo(f'Issue {issue.id}: {issue.title}')
        s = assess_issue_completeness(issue=issue, graph=graph)
        d = s.summary_dict()
        yield d


def assess_issue_completeness(*, issue: Issue, graph: BELGraph) -> CompletenessSummary:
    """Summarize the novelty of the PMID referenced by the issue."""
    pmid = _get_pmid(issue.description)
    ids = ('pmid', pmid)
    return assess_completeness(ids, graph)


def _get_pmid(description: str) -> Optional[str]:
    for line in description.split('\n'):
        line = line.strip()
        if line.startswith(_prefix):
            line: str = line[len(_prefix):]
            line = line[:line.index(']')]
            return line


if __name__ == '__main__':
    main()
