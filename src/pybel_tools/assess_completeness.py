# -*- coding: utf-8 -*-

"""Assess the completeness of a graph.

Run on CONIB with ``python -m pybel_tools.assess_completeness [PMID]``.
"""

import logging
import math
from dataclasses import dataclass
from typing import Iterable, List, Optional, Set, TextIO, Tuple, Union

import click
from indra.sources.indra_db_rest.api import get_statements_for_paper
from indra.sources.indra_db_rest.util import logger as indra_logger

from pybel import BELGraph, BaseEntity, from_indra_statements
from pybel.cli import _from_pickle_callback

__all__ = [
    'assess_completeness',
    'CompletenessSummary',
]


@dataclass
class CompletenessSummary:
    """Information about the assessment of the completeness of a reference graph versus a set of documents."""

    documents: List[Tuple[str, str]]
    reference_nodes: Set[BaseEntity]
    indra_nodes: Set[BaseEntity]

    @property
    def novel_nodes(self) -> Set[BaseEntity]:  # noqa: D401
        """The nodes INDRA found that weren't in the reference."""
        return self.indra_nodes - self.reference_nodes

    def summary_dict(self):
        """Summarize the findings in a tuple."""
        curies = [
            f'{x}:{y}'
            for x, y in self.documents
        ]
        indra_novel = len(self.novel_nodes) / len(self.indra_nodes)
        reference_novel = len(self.novel_nodes) / len(self.reference_nodes)

        # normalize by the number of nodes and the number of documents
        score_novel = len(self.novel_nodes) / math.sqrt(len(self.indra_nodes) * len(self.reference_nodes)) / len(
            self.documents)

        return dict(
            curies=curies,
            indra_count=len(self.indra_nodes),
            indra_novel=indra_novel,
            reference_count=len(self.reference_nodes),
            reference_novel=reference_novel,
            nodel_count=len(self.novel_nodes),
            score_novel=score_novel,
        )

    def summary_str(self) -> str:
        """Summarize the findings in a string."""
        d = self.summary_dict()
        return f"""
Novelty check for {', '.join(d['curies'])}:
    Reference had: {d['reference_count']} nodes
    INDRA found: {d['indra_count']} nodes ({d['novel_count']} new, {d['indra_novel']:.2%})
    Curation novelty: {d['reference_novel']:.2%}
    Score: {d['score_novel']:.2%}
"""

    def summarize(self, file: Optional[TextIO] = None) -> None:
        """Print the summary of the findings."""
        print(self.summary_str(), file=file)


def assess_completeness(
    ids: Union[str, Tuple[str, str], List[Tuple[str, str]]],
    reference_nodes: Union[Iterable[BaseEntity], BELGraph],
    *,
    verbose_indra_logger: bool = False,
) -> Optional[CompletenessSummary]:
    """Check INDRA if the given document has new interesting content compared to the graph.

    :param ids: A CURIE (e.g., pmid:30606258, pmc:PMC6318896),
     a pair of database/identifier (e.g. `('pmid', '30606258')`) or a list of pairs.
    :param reference_nodes: A set of nodes or a BEL graph (which is a set of nodes)
    :param verbose_indra_logger: Should the INDRA logger show any output below the WARNING level?
    """
    if not verbose_indra_logger:
        indra_logger.setLevel(logging.WARNING)

    if isinstance(ids, str):
        _ids = ids.split(':')
        if len(_ids) != 2:
            raise ValueError(f'String should be given as CURIE: {ids}')
        ids = [_ids]
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
        return

    return CompletenessSummary(
        documents=ids,
        reference_nodes=set(reference_nodes),
        indra_nodes=indra_nodes,
    )


@click.command()
@click.argument('pmids', nargs=-1)
@click.option(
    '-g', '--graph',
    metavar='path',
    type=click.File('rb'),
    callback=_from_pickle_callback,
    help='Path to BEL file. Loads CONIB as an example if none given',
)
@click.option('-v', '--verbose', is_flag=True)
def main(pmids: str, graph: Optional[BELGraph], verbose: bool) -> None:
    """Check a BEL graph for added value of a given article.

    Example: 30606258 for paper entitled "A pathogenic tau fragment compromises microtubules,
    disrupts insulin signaling and induces the unfolded protein response."
    """
    if graph is None:
        import hbp_knowledge
        graph = hbp_knowledge.get_graph()

    pmids = [('pmid', pmid) for pmid in pmids]

    s = assess_completeness(pmids, graph, verbose_indra_logger=verbose)
    s.summarize()


if __name__ == '__main__':
    main()
