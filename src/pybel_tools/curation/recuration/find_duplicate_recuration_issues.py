# -*- coding: utf-8 -*-

"""A script that finds duplicate re-curation issues."""

from collections import defaultdict
from typing import List, Mapping

import click
from easy_config.contrib.click import args_from_config
from gitlab.v4.objects import Issue, Project
from pybel_git.gitlab import GitlabConfig

from .utils import RECURATION_ISSUE_NAME_PREFIX, RECURATION_LABEL


def get_recuration_issues(project) -> Mapping[str, List[Issue]]:
    issues = defaultdict(list)
    for recuration_issue in project.issues.list(all=True, labels=[RECURATION_LABEL]):
        if recuration_issue.title.lower().startswith(RECURATION_ISSUE_NAME_PREFIX):
            parts = recuration_issue.title[len(RECURATION_ISSUE_NAME_PREFIX):].split()
            key = parts[0].strip().lower()
            issues[key].append(recuration_issue)
    return issues


@click.command()
@args_from_config(GitlabConfig)
def main(project_id: int, url: str, token: str):
    gitlab_config = GitlabConfig.load(  # noqa: S106
        project_id=project_id,
        url=url,
        token=token,
    )
    project = gitlab_config.get_project()
    _do_it(project)


def _do_it(project: Project) -> None:
    for title, issues in get_recuration_issues(project).items():
        if len(issues) > 1:
            print(title, issues)
            for issue in issues:
                print(f'https://gitlab.scai.fraunhofer.de/charles.hoyt/hbp/issues/{issue.get_id()}')


if __name__ == '__main__':
    main()
