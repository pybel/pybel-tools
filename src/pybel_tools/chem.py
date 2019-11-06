# -*- coding: utf-8 -*-

"""Chemistry tools for BEL."""

import itertools as itt
from typing import Iterable, Tuple

from rdkit import DataStructs
from rdkit.Chem import MACCSkeys, MolFromInchi
from tqdm import tqdm

from pybel import BELGraph
from pybel.constants import ANNOTATIONS, IDENTIFIER, NAMESPACE
from pybel.dsl import BaseEntity
from pybel.struct import invert_node_predicate, node_predicate, remove_filtered_nodes

__all__ = [
    'add_similarity_edges',
    'remove_non_inchi',
]


def iter_inchi_nodes(graph: BELGraph) -> Iterable[Tuple[tuple, dict, str]]:
    """Iterate over node tuple, node data, and InChI string triples."""
    for node, data in graph.nodes(data=True):
        if node_is_inchi(data):
            yield node, data, data.get(IDENTIFIER)


@node_predicate
def node_is_inchi(node: BaseEntity):
    return 'inchi' == node.get(NAMESPACE) and IDENTIFIER in node


def remove_non_inchi(graph: BELGraph):
    """Remove all non-inchi nodes."""
    remove_filtered_nodes(graph, invert_node_predicate(node_is_inchi))


def add_similarity_edges(graph: BELGraph, cutoff: float = 0.8) -> None:
    """Enrich a BEL graph with edges between chemicals that have InChI.

    :param graph: A BEL graph
    :param cutoff: The cutoff for similarity
    """
    inchi_to_node_tuple = {
        inchi: node_tuple
        for node_tuple, _, inchi in iter_inchi_nodes(graph)
    }

    mols = {
        inchi: MolFromInchi(inchi)
        for inchi in inchi_to_node_tuple
    }

    fps = {
        inchi: MACCSkeys.GenMACCSKeys(mol)
        for inchi, mol in tqdm(mols.items(), desc='calculating MACCS keys')
        if mol is not None
    }

    n_combinations = (len(fps) * (len(fps) - 1) / 2)
    _sim_iter = (
        (x, y, DataStructs.FingerprintSimilarity(fps[x], fps[y]))
        # can also use FingerprintMols
        for x, y in tqdm(itt.combinations(fps, 2), total=n_combinations, desc='calculating similarity')
    )

    for x, y, sim in _sim_iter:
        if sim < cutoff:
            continue

        source, target = inchi_to_node_tuple[x], inchi_to_node_tuple[y]

        key = graph.add_unqualified_edge(source, target, relation='similar')
        graph[source][target][key][ANNOTATIONS] = {
            'similarity': sim
        }
