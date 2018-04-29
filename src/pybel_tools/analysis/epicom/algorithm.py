# -*- coding: utf-8 -*-

"""This algorithm has multiple steps

1. Select NeuroMMSig networks for AD, PD, and epilepsy
2. Select drugs from DrugBank, and their targets
3. Run NeuroMMSig algorithm on target list for each network and each mechanism
4. Store in database
"""

import itertools as itt
import logging
import os

import bio2bel_drugbank
from bio2bel_drugbank.constants import DATA_DIR as DRUGBANK_DATA_DIR
from tqdm import tqdm

from pybel.dsl import gene as gene_dsl
from pybel_tools.analysis.neurommsig import get_neurommsig_score, neurommsig_graph_preprocessor
from pybel_tools.grouping import get_subgraphs_by_annotation
from pybel_tools.summary import get_annotation_values

log = logging.getLogger(__name__)

_dti_cache_path = os.path.join(DRUGBANK_DATA_DIR, 'drug_to_gene_symbols.json')


def _get_drug_target_interactions():
    """Gets a mapping from drugs to their list of gene

    :rtype: dict[str,list[str]]
    """
    drugbank_manager = bio2bel_drugbank.Manager()
    if not drugbank_manager.is_populated():
        drugbank_manager.populate()

    return drugbank_manager.get_drug_to_hgnc_symbols()


def _preprocess_dtis(dtis):
    return {
        drug: [gene_dsl(namespace='HGNC', name=target).as_tuple() for target in targets]
        for drug, targets in dtis.items()
    }


def epicom_on_graph(graph, dtis, preprocess=True):
    """
    :param pybel.BELGraph graph:
    :param dict[str,list[tuple]] dtis:
    :param bool preprocess: If true, preprocess the graph with :func:`neurommsig_graph_preprocessor`.
    :rtype: iter[tuple[str,str,float]]
    """
    if preprocess:
        log.info('preprocessing %s', graph)
        graph = neurommsig_graph_preprocessor(graph)

    log.info('stratifying %s', graph)
    subgraphs = get_subgraphs_by_annotation(graph, annotation='Subgraph', keep_undefined=False)

    log.info('running subgraphs x drugs for %s', graph)
    it = itt.product(sorted(subgraphs), sorted(dtis))
    it = tqdm(it, total=len(subgraphs) * len(dtis), desc='Calculating scores')

    def get_score(s, d):
        """Gets the score

        :param str s: name of the subgraph
        :param str d: name of the drug
        :rtype: Optional[float]
        """
        return get_neurommsig_score(subgraphs[s], dtis[d])

    for subgraph_name, drug in it:
        score = get_score(subgraph_name, drug)

        if score is None or score == 0.0:
            continue

        yield drug, subgraph_name, score


def _multi_run_helper(graphs):
    """
    :param iter[pybel.BELGraph] graphs: A BEL Graph
    :param dict[str,list[str]] dtis: A dictionary of drugs mapping to their targets
    :rtype: iter[tuple[str,str,str,float]]
    """
    dtis = _preprocess_dtis(_get_drug_target_interactions())

    for graph in graphs:
        for drug, subgraph_name, score in epicom_on_graph(graph, dtis):
            yield graph.name, subgraph_name, drug, score


def _multi_run_helper_file_wrapper(graphs, file):
    for row in _multi_run_helper(graphs):
        print(*row, sep='\t', file=file)


def multi_run_epicom(graphs, path):
    """Run EpiCom analysis  on many graphs

    :param iter[pybel.BELGraph] graphs:
    :param str or file path: output file path
    """
    if isinstance(path, str):
        with open(path, 'w') as file:
            _multi_run_helper_file_wrapper(graphs, file)

    else:
        _multi_run_helper_file_wrapper(graphs, path)


def run_epicom(graph, directory):
    """
    :param pybel.BELGraph graph: A BEL Graph
    :param str directory: The directory in which the algorithm is run
    :rtype: iter[tuple[str,str,str,float]]
    """
    os.makedirs(directory, exist_ok=True)

    dtis = _preprocess_dtis(_get_drug_target_interactions())

    drugs = list(dtis)
    subgraphs = get_annotation_values(graph, annotation='Subgraph')

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
        for drug, subgraph_name, score in epicom_on_graph(graph, dtis):
            print(
                subgraph_name_to_id[subgraph_name],
                drug_name_to_id[drug],
                score,
                sep='\t',
                file=file
            )
