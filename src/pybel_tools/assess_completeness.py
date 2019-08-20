# -*- coding: utf-8 -*-

"""Assess the completeness of a graph.

Run on CONIB with ``python -m pybel_tools.assess_completeness [PMID]``.
"""

import logging

import click
from indra.sources.indra_db_rest.api import get_statements_for_paper
from indra.sources.indra_db_rest.util import logger as indra_logger

from pybel import BELGraph, from_indra_statements

__all__ = [
    'assess_completeness',
]


def assess_completeness(ids, graph: BELGraph, verbose_indra_logger: bool = False):
    """Check INDRA if the given document has new interesting content compared to the graph.

    :param ids: A CURIE (e.g., pmid:30606258, pmc:PMC6318896),
     a pair of database/identifier (e.g. `('pmid', '30606258')`) or a list of pairs.
    :param graph: A BEL graph
    """
    if not verbose_indra_logger:
        indra_logger.setLevel(logging.WARNING)

    if isinstance(ids, str):
        ids = [ids.split(':')]
    elif isinstance(ids, tuple):
        ids = [ids]

    # Normalize PMC database name as well as stringify all identifiers
    ids = [
        ('pmcid' if db == 'pmc' else db, str(db_id))
        for db, db_id in ids
    ]

    stmts = get_statements_for_paper(ids=ids)
    indra_graph = from_indra_statements(stmts)

    indra_nodes = set(indra_graph)
    if not indra_nodes:
        print(f'INDRA did not return any results for {ids}')
        return False

    query_nodes = set(graph)
    new_nodes = indra_nodes - query_nodes

    print(f"""
    Graph had:
        Total nodes: {len(query_nodes)}

    INDRA found (in {ids}):
        Total nodes: {len(indra_nodes)}
        New nodes: {len(new_nodes)} ({len(new_nodes) / len(indra_nodes):.2%})
    """)
    return True


@click.command()
@click.argument('pmid')
@click.option('-v', '--verbose', is_flag=True)
def main(pmid: str, verbose: bool) -> None:
    """Check CONIB for added value of a given article.

    Example: 30606258 for paper entitled "A pathogenic tau fragment compromises microtubules,
    disrupts insulin signaling and induces the unfolded protein response."
    """
    import hbp_knowledge
    graph = hbp_knowledge.get_graph()

    assess_completeness(('pmid', pmid), graph, verbose_indra_logger=verbose)


if __name__ == '__main__':
    main()
