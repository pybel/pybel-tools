# -*- coding: utf-8 -*-

"""An implementation of a drug-target-based mechanism enrichment strategy."""

import itertools as itt
import logging
import os
from typing import Iterable, List, Mapping, Optional, TextIO, Tuple, Union

from tqdm import tqdm

from pybel import BELGraph
from pybel.dsl import Gene
from pybel.struct.grouping import get_subgraphs_by_annotation
from pybel.struct.summary import get_annotation_values
from ..neurommsig import get_neurommsig_score, neurommsig_graph_preprocessor

__all__ = [
    'run_epicom',
    'get_drug_scores',
]

logger = logging.getLogger(__name__)


def run_epicom(graph: BELGraph, directory: str, annotation: str = 'Subgraph') -> None:
    """Run EpiCom reloaded on the given graph stratifying with the given annotation.

    :param graph: A BEL graph
    :param directory: The location to output the results
    :param annotation: The annotation to use to stratify the graph into sub-graphs.
    """
    os.makedirs(directory, exist_ok=True)

    dtis = _preprocess_dtis(_get_drug_target_interactions())

    drugs = list(dtis)
    subgraphs = get_annotation_values(graph, annotation=annotation)

    subgraph_name_to_id = {name: i for i, name in enumerate(sorted(subgraphs))}
    drug_name_to_id = {name: i for i, name in enumerate(sorted(drugs))}

    with open(os.path.join(directory, 'subgraphs.tsv'), 'w') as file:
        print('id', 'name', sep='\t', file=file)
        for i, name in enumerate(sorted(subgraphs)):
            print(i, name, sep='\t', file=file)

    with open(os.path.join(directory, 'drugs.tsv'), 'w') as file:
        print('id', 'name', sep='\t', file=file)
        for i, name in enumerate(sorted(dtis)):
            print(i, name, sep='\t', file=file)

    with open(os.path.join(directory, 'scores.tsv'), 'w') as file:
        print('subgraph', 'drug', 'score', sep='\t', file=file)
        for drug, subgraph_name, score in get_drug_scores(graph, dtis):
            print(
                subgraph_name_to_id[subgraph_name],
                drug_name_to_id[drug],
                score,
                sep='\t',
                file=file
            )


def get_drug_scores(
        graph: BELGraph,
        dtis: Mapping[str, List[Gene]],
        preprocess_graph: bool = True,
) -> Iterable[Tuple[str, str, float]]:
    """Get drug scores for the given graph."""
    if preprocess_graph:
        logger.info('preprocessing %s', graph)
        graph = neurommsig_graph_preprocessor(graph)

    logger.info('stratifying %s', graph)
    subgraphs = get_subgraphs_by_annotation(graph, annotation='Subgraph', sentinel='UNDEFINED')

    logger.info('running subgraphs x drugs for %s', graph)
    it = itt.product(sorted(subgraphs), sorted(dtis))
    it = tqdm(it, total=len(subgraphs) * len(dtis), desc='Calculating scores')

    def get_score(s: str, d: str) -> Optional[float]:
        """Get the score for the given subgraph and drug."""
        return get_neurommsig_score(subgraphs[s], dtis[d])

    for subgraph_name, drug in it:
        score = get_score(subgraph_name, drug)

        if score is None or score == 0.0:
            continue

        yield drug, subgraph_name, score


def _get_drug_target_interactions(manager: Optional['bio2bel_drugbank.manager'] = None) -> Mapping[str, List[str]]:
    """Get a mapping from drugs to their list of gene."""
    if manager is None:
        import bio2bel_drugbank
        manager = bio2bel_drugbank.Manager()

    if not manager.is_populated():
        manager.populate()

    return manager.get_drug_to_hgnc_symbols()


def _preprocess_dtis(dtis: Mapping[str, List[str]]) -> Mapping[str, List[Gene]]:
    return {
        drug: [Gene(namespace='HGNC', name=target) for target in targets]
        for drug, targets in dtis.items()
    }


def _multi_run_helper(graphs: Iterable[BELGraph]) -> Iterable[Tuple[str, str, str, float]]:
    dtis = _preprocess_dtis(_get_drug_target_interactions())

    for graph in graphs:
        for drug, subgraph_name, score in get_drug_scores(graph, dtis):
            yield graph.name, subgraph_name, drug, score


def _multi_run_helper_file_wrapper(graphs: Iterable[BELGraph], file: Optional[TextIO] = None) -> None:
    for row in _multi_run_helper(graphs):
        print(*row, sep='\t', file=file)


def multi_run_epicom(graphs: Iterable[BELGraph], path: Union[None, str, TextIO]) -> None:
    """Run EpiCom analysis on many graphs."""
    if isinstance(path, str):
        with open(path, 'w') as file:
            _multi_run_helper_file_wrapper(graphs, file)
    else:
        _multi_run_helper_file_wrapper(graphs, path)
