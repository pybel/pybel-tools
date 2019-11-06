# -*- coding: utf-8 -*-

"""Utility functions for re-curation tools."""

from typing import Mapping

from gitlab.v4.objects import Issue

CURATION_LABEL = 'Curation'
RECURATION_LABEL = 'Re-curation'
RECURATION_ISSUE_NAME_PREFIX = 're-curate '.lower()
CURATION_ISSUE_NAME_PREFIX = 'curate '.lower()
RECURATION_ISSUE_MISSING_DESCRIPTION = 'Missing original issue.'


def get_recuration_issues(project) -> Mapping[str, Issue]:
    return _get_issue_prefixed(project, RECURATION_LABEL, RECURATION_ISSUE_NAME_PREFIX)


def get_curation_issues(project, all_issues: bool = True) -> Mapping[str, Issue]:
    return _get_issue_prefixed(project, CURATION_LABEL, CURATION_ISSUE_NAME_PREFIX, all_issues=all_issues)


def _get_issue_prefixed(project, label, prefix, all_issues: bool = True):
    issues = {}
    for issue in project.issues.list(all=all_issues, labels=[label]):
        if issue.title.lower().startswith(prefix):
            title = issue.title[len(prefix):]
            parts = title.split()
            key = parts[0].strip().strip(':').lower()
            issues[key] = issue
    return issues
