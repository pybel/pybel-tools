# -*- coding: utf-8 -*-

"""Curation tools for the Gene Ontology (GO).

Run with `python -m hbp.curation.planning.pathways`.
"""

from typing import Iterable, Optional

import click

import bio2bel_kegg
import bio2bel_reactome
import bio2bel_wikipathways
from compath_utils import CompathManager
from pybel.cli import connection_option


def get_managers(connection: Optional[str] = None) -> Iterable[CompathManager]:
    wikipathways_manager = bio2bel_wikipathways.Manager(connection=connection)
    if not wikipathways_manager.is_populated():
        click.echo('WikiPathways is not populated')
    else:
        yield wikipathways_manager

    reactome_manager = bio2bel_reactome.Manager(connection=connection)
    if not reactome_manager.is_populated():
        click.echo('Reactome is not populated')
    else:
        yield reactome_manager

    kegg_manager = bio2bel_kegg.Manager(connection=connection)
    if not kegg_manager.is_populated():
        click.echo('KEGG is not populated')
    else:
        yield kegg_manager


@click.command()
@click.argument('name')
@click.option('--namespace', type=click.Choice('hgnc.symbol'), default='hgnc.symbol')
@connection_option
def main(name: str, namespace: str, connection: Optional[str]):
    for manager in get_managers(connection):
        if namespace == 'hgnc.symbol':
            protein = manager.get_protein_by_hgnc_symbol(name)
        else:
            raise ValueError(f'{namespace} is not yet supported')

        if protein is None:
            click.echo(f'No pathways in {manager.module_name}')
        else:
            for pathway in protein.pathways:
                pathway_id = getattr(pathway, f'{manager.module_name}_id')
                click.echo(f'{manager.module_name}:{pathway_id} ! {pathway}')


if __name__ == '__main__':
    main()
