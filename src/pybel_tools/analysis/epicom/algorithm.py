# -*- coding: utf-8 -*-

"""This algorithm has multiple steps

1. Select NeuroMMSig networks for AD, PD, and epilepsy
2. Select drugs from DrugBank, and their targets
3. Run NeuroMMSig algorithm on target list for each network and each mechanism
4. Store in database
"""

import itertools as itt
import json
import logging
import os
from collections import defaultdict

from bio2bel_drugbank import Manager as DrugBankManager
from bio2bel_drugbank.constants import DATA_DIR as DRUGBANK_DATA_DIR
from tqdm import tqdm

from bio2bel_hgnc import Manager as HGNCManager
from pybel.dsl import gene as gene_dsl
from pybel_tools.analysis.neurommsig import get_neurommsig_score, neurommsig_graph_preprocessor
from pybel_tools.grouping import get_subgraphs_by_annotation

log = logging.getLogger(__name__)

_dti_cache_path = os.path.join(DRUGBANK_DATA_DIR, 'drug_genes.json')


def _get_drug_target_interactions():
    """Gets a mapping from drugs to their list of gene

    :rtype: dict[str,list[str]]
    """

    if os.path.exists(_dti_cache_path):
        log.info('loading cached DTIs')
        with open(_dti_cache_path) as file:
            return json.load(file)

    hgnc_manager = HGNCManager()
    if not hgnc_manager.is_populated():
        hgnc_manager.populate()

    hgnc_id_symbol_mapping = hgnc_manager.build_hgnc_id_symbol_mapping()

    drugbank_manager = DrugBankManager()
    if not drugbank_manager.is_populated():
        drugbank_manager.populate()

    rv = defaultdict(list)

    for dpi in tqdm(drugbank_manager.list_drug_protein_interactions(),
                    total=drugbank_manager.count_drug_protein_interactions(),
                    desc='getting DTIs'):

        if dpi.protein.hgnc_id is None:
            continue

        hgnc_id = dpi.protein.hgnc_id[len('HGNC:'):] if dpi.protein.hgnc_id.startswith('HGNC:') else dpi.protein.hgnc_id

        hgnc_symbol = hgnc_id_symbol_mapping.get(hgnc_id)

        if hgnc_symbol is None:
            log.warning('could not map HGNC identifier: %s', hgnc_id)

        rv[dpi.drug.name].append(hgnc_symbol)

    with open(_dti_cache_path, 'w') as file:
        log.info('dumping cached DTIs')
        json.dump(rv, file)

    return rv


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
    subgraphs = get_subgraphs_by_annotation(graph, annotation='Subgraph')

    log.info('running subgraphs x drugs for %s', graph)
    it = itt.product(subgraphs.items(), dtis.items())
    it = tqdm(it, total=len(subgraphs) * len(dtis))

    for (subgraph_name, subgraph), (drug, genes) in it:
        score = get_neurommsig_score(subgraph, genes)

        if score is None:
            continue

        yield drug, subgraph_name, score


def _run_helper(graphs):
    """
    :param iter[pybel.BELGraph] graphs: A list of BEL grahs
    :param dict[str,list[str]] dtis: A dictionary of drugs mapping to their targets
    :rtype: iter[tuple[str,str,str,float]]
    """
    dtis = _preprocess_dtis(_get_drug_target_interactions())

    for graph in graphs:
        for drug, subgraph_name, score in epicom_on_graph(graph, dtis):
            yield graph.name, subgraph_name, drug, score


def _run_helper_file_wrapper(graphs, file):
    for row in _run_helper(graphs):
        print(*row, sep='\t', file=file)


def run_epicom(graphs, path):
    """Run EpiCom analysis  on many graphs

    :param iter[pybel.BELGraph] graphs:
    :param str or file path: output file path
    """
    if isinstance(path, str):
        with open(path, 'w') as file:
            _run_helper_file_wrapper(graphs, file)

    else:
        _run_helper_file_wrapper(graphs, path)
