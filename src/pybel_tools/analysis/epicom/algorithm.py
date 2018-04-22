# -*- coding: utf-8 -*-

"""This algorithm has multiple steps

1. Select NeuroMMSig networks for AD, PD, and epilepsy
2. Select drugs from DrugBank, and their targets
3. Run NeuroMMSig algorithm on target list for each network and each mechanism
4. Store in database
"""

from collections import defaultdict

from bio2bel_drugbank import Manager as DrugBankManager

from pybel_tools.analysis.neurommsig import get_neurommsig_scores_stratified


def _get_drug_target_interactions():
    """Gets a mapping from drugs to their list of targets

    :rtype: dict[pybel.dsl.abundance,list[pybel.dsl.protein]]
    """
    manager = DrugBankManager()
    rv = defaultdict(list)

    for dpi in manager.list_drug_protein_interactions():
        rv[dpi.drug.name].append(dpi.protein.hgnc_id)

    return rv


def _run_graph(graph, dtis):
    rv = {}

    for drug, targets in dtis.items():
        scores = get_neurommsig_scores_stratified(graph, targets)
        for subgraph_name, score in scores.items():
            rv[drug, subgraph_name] = score

    return rv


def _run_helper(graphs, dtis):
    """
    :param iter[pybel.BELGraph] graphs: A list of BEL grahs
    :param dict[pybel.dsl.abundance,list[pybel.dsl.protein]] dtis: A dictionary of drugs mapping to their targets
    """
    rv = {}

    for graph in graphs:
        scores = _run_graph(graph, dtis)
        for (drug, subgraph_name), score in scores.items():
            rv[graph.name, subgraph_name, drug] = score

    return rv


def run_epicom(graphs):
    dtis = _get_drug_target_interactions()
    return _run_helper(graphs, dtis)
