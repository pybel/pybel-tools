# -*- coding: utf-8 -*-

"""Curation tools for PubMed.

Run with `python -m hbp.curation.planning.pubmed`.
"""

import sys
from typing import List

import click
from easy_config.contrib.click import args_from_config
from pybel_git.gitlab import GitlabConfig

from ..utils import make_issues_from_pmids, min_year_option


@click.command()
@args_from_config(GitlabConfig)
@click.option('-f', '--file', default=sys.stdin, type=click.File())
@min_year_option
@click.option('--make-issues', is_flag=True, help='Create issues on GitLab HBP repository')
@click.option('--allow-closed', is_flag=True, help='Allow publications that are not on PMC')
@click.option('-l', '--label', multiple=True)
def main(project_id: int, url: str, token: str, file, min_year: int, make_issues: bool, allow_closed: bool,
         label: List[str]):
    """Get a list of documents by their PubMed identifiers."""
    gitlab_config = GitlabConfig.load(  # noqa: S106
        project_id=project_id,
        url=url,
        token=token,
    )
    project = gitlab_config.get_project()

    pmids = list(sorted({
        line.strip()
        for line in file
    }))

    make_issues_from_pmids(
        project,
        pmids,
        min_year=min_year,
        allow_closed=allow_closed,
        make_issues=make_issues,
        labels=label,
    )


if __name__ == '__main__':
    main()
