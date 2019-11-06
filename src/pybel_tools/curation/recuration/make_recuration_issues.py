# -*- coding: utf-8 -*-

"""A one-off script for making issues in GitLab for BEL graphs without re-curation.

Run with: ``python -m hbp.curation.recuration.make_recuration_issues``.

If you just want to test what would happen, run with ``--dry``.
"""

import os

import click
import hbp_knowledge
from easy_config.contrib.click import args_from_config
from pybel import BELGraph
from pybel.constants import ANNOTATIONS, RELATION, UNQUALIFIED_EDGES
from pybel_git.gitlab import GitlabConfig

from .utils import (
    CURATION_LABEL, RECURATION_ISSUE_MISSING_DESCRIPTION,
    RECURATION_ISSUE_NAME_PREFIX, RECURATION_LABEL, get_curation_issues, get_recuration_issues,
)


@click.command()
@args_from_config(GitlabConfig)
@click.option('--dry', is_flag=True)
@click.option('-v', '--verbose', is_flag=True)
def main(project_id: int, url: str, token: str, dry: bool, verbose: bool):
    """Reorganize re-curation issues."""
    gitlab_config = GitlabConfig.load(  # noqa: S106
        project_id=project_id,
        url=url,
        token=token,
    )
    project = gitlab_config.get_project()

    click.echo('Getting re-curation issues')
    recuration_issues = get_recuration_issues(project)
    click.echo(f'Retrieved {len(recuration_issues)} re-curation issues labeled "{RECURATION_LABEL}"')

    click.echo('Getting curation issues')
    original_issues = get_curation_issues(project)
    click.echo(f'Retrieved {len(original_issues)} curation issues labeled "{CURATION_LABEL}"')

    overlap = set(recuration_issues).intersection(original_issues)
    click.echo(f'Overlap: {len(overlap)}')

    for i in recuration_issues:
        if i not in original_issues:
            print(i)
    for i in original_issues:
        if i not in recuration_issues:
            print(i)

    click.echo('Getting graphs')
    graphs = hbp_knowledge.get_graphs()
    click.echo(f'Got {len(graphs)} graphs')

    new_count, closed_count, done_count, pending_count = 0, 0, 0, 0
    for path, graph in graphs.items():
        recurated = graph_has_recuration(graph)
        basename = os.path.basename(path)[:-len('.bel')]
        recuration_issue = recuration_issues.get(basename)

        original_issue = original_issues.get(basename)
        ss = '\n' + '\n- '.join(graph.summary_str().split('\n'))
        if original_issue:
            description = f'Original issue: #{original_issue.get_id()}\n{ss}'
        else:
            description = RECURATION_ISSUE_MISSING_DESCRIPTION + f'\n{ss}'

        if recuration_issue is None and recurated and verbose:
            click.echo(f'✅ {basename} already re-curated')

        elif recuration_issue is None and not recurated:
            title = make_issue_title(basename)
            click.echo(f'❌ {basename} creating issue: {title}')
            new_count += 1
            if not dry:
                recuration_issue = project.issues.create({
                    'title': title,
                    'description': description,
                })
                recuration_issue.labels.append(RECURATION_LABEL)
                recuration_issue.save()

        elif recuration_issue and recurated:
            if recuration_issue.state == 'opened':
                click.echo(f'🚪 {basename} closing issue {recuration_issue.title}')
                closed_count += 1
                if not dry:
                    recuration_issue.state_event = 'close'
                    recuration_issue.save()
            else:
                if verbose:
                    click.echo(f'✅ {basename} already done')
                done_count += 1

        elif recuration_issue and not recurated:
            click.echo(f'❌ {basename} not yet re-curated (#{recuration_issue.get_id()})')
            pending_count += 1

        if recuration_issue is not None and graph.description != description:
            click.echo(f'{basename} updating description')
            if not dry:
                recuration_issue.description = description
                recuration_issue.save()

    click.echo(f'''
Summary
===========
    New {new_count}
 Closed {closed_count}
Pending {pending_count}
   Done {done_count}
''')


def make_issue_title(name: str) -> str:
    """Format the graph's name into an issue name."""
    return f'{RECURATION_ISSUE_NAME_PREFIX}{name}'


def graph_has_recuration(graph: BELGraph) -> bool:
    """Check that all edges have re-curation."""
    return all(
        ANNOTATIONS in data and 'Confidence' in data[ANNOTATIONS]
        for _, _, data in graph.edges(data=True)
        if data[RELATION] not in UNQUALIFIED_EDGES
    )


if __name__ == '__main__':
    main()
