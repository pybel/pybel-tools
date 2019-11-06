# -*- coding: utf-8 -*-

"""Curation tools for the Gene Ontology (GO).

Run with `python -m hbp.curation.planning.go`.
"""

from typing import List

import click
import requests

from ..utils import make_issues_from_pmids, min_year_option

url = 'http://golr-aux.geneontology.io/solr/select'

BASE_PARAMS = {
    'defType': ['edismax'],
    'qt': ['standard'],
    'indent': ['on'],
    'wt': ['csv'],
    'rows': ['100000'],
    'start': ['0'], 'fl': ['reference'],
    'facet': ['true'],
    'facet.mincount': ['1'],
    'facet.sort': ['count'],
    'json.nl': ['arrarr'],
    'facet.limit': ['25'],
    'hl': ['true'],
    'hl.simple.pre': ['<em class="hilite">'],
    'hl.snippets': ['1000'],
    'csv.separator': ['\t'],
    'csv.header': ['false'],
    'csv.mv.separator': ['|'],
    'fq': ['document_category:"annotation"'],  # add bioentity here too
    'facet.field': ['aspect', 'taxon_subset_closure_label', 'type', 'evidence_subset_closure_label',
                    'regulates_closure_label', 'annotation_class_label', 'qualifier',
                    'annotation_extension_class_closure_label', 'assigned_by', 'panther_family_label'],
    'q': ['*:*'],
}


def get_pmids_from_go_annotations_by_uniprot_id(uniprot_id: str) -> List[str]:
    """Get the PubMed identifiers used in GO annotations for the given protein."""
    params = BASE_PARAMS.copy()
    params['fq'].append(f'bioentity:"UniProtKB:{uniprot_id}"')
    r = requests.get(url, params)
    lines = (
        line.strip()
        for line in r.text.splitlines()
    )
    return list(sorted({
        line.split(':')[1]
        for line in lines
        if line and line.lower().startswith('pmid')
    }))


@click.command()
@click.argument('uniprot_id')
@click.option('--namespace', type=click.Choice(['uniprot']), default='uniprot')
@min_year_option
@click.option('--make-issues', is_flag=True, help='Create issues on GitLab HBP repository')
@click.option('--allow-closed', is_flag=True, help='Allow publications that are not on PMC')
@click.option('-l', '--label', multiple=True)
def main(uniprot_id: str, namespace: str, min_year: int, make_issues: bool, allow_closed: bool, label: List[str]):
    """Get a list of documents for the given UniProt identifier.

    Example: Q13148.
    """
    if namespace == 'uniprot':
        pmids = get_pmids_from_go_annotations_by_uniprot_id(uniprot_id)
    else:
        raise ValueError(f'{namespace} is not yet supported')

    make_issues_from_pmids(
        pmids,
        min_year=min_year,
        allow_closed=allow_closed,
        make_issues=make_issues,
        labels=label,
    )


if __name__ == '__main__':
    main()
