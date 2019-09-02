# -*- coding: utf-8 -*-

"""Chemistry utilities for BEL that require RDKit."""

import itertools as itt
from typing import Iterable, Mapping, Tuple

import requests
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys
from tqdm import tqdm

from pybel import Abundance, BELGraph, BaseAbundance

__all__ = [
    'add_chemical_similarities',
]


def add_chemical_similarities(graph: BELGraph, cutoff: float = 0.7) -> None:
    """Add chemical similarities to the graph."""
    node_to_smiles = get_node_to_smiles(graph)
    similarities = get_similarity(node_to_smiles)
    for (source, target), similarity in tqdm(similarities.items(), desc='Updating graph'):
        if cutoff < similarity:
            graph.add_unqualified_edge(source, target, 'chemical_similarity')


def get_node_to_smiles(graph: BELGraph) -> Mapping[BaseAbundance, str]:
    """Get a mapping from nodes to SMILES strings."""
    return dict(_help_iter_node_smiles(graph))


def _help_iter_node_smiles(graph: BELGraph) -> Iterable[Tuple[Abundance, str]]:
    for node in graph:
        if not isinstance(node, Abundance):
            continue
        if node.namespace.lower() == 'pubchem.compound':
            yield node, cid_to_smiles(node.identifier)
        elif node.namespace.lower() == 'smiles':
            yield node, node.identifier


def cid_to_smiles(pubchem_id: str) -> str:
    """Get the SMILES for chemicals in PubChem database.

    :param pubchem_id: PubChem compound identifier
    :return: SMILES
    """
    url = f"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_id}/property/canonicalSMILES/TXT"
    response = requests.get(url)
    return response.text


def get_similarity(node_to_smiles: Mapping[BaseAbundance, str]):
    """Get the similarities between all pair combinations of chemicals in the list.

    :param node_to_smiles: a dictionary with nodes as keys and smiles as values
    :return: a dictionary with the pair chemicals as key and similarity calculation as value
    """
    fps = get_fingerprints(node_to_smiles)
    return {
        (source, target): DataStructs.FingerprintSimilarity(mol_1, mol_2)
        for (source, mol_1), (target, mol_2) in
        tqdm(itt.combinations(fps.items(), 2), desc='Calculating Similarities')
    }


def get_fingerprints(node_to_smiles: Mapping[BaseAbundance, str]):
    """Create a dictionary containing the fingerprints for every chemical.

    :param node_to_smiles: a dictionary with nodes as keys and smiles as values
    :return: a dictionary with nodes as keys and the MACCSkeys fingerprints
    """
    node_to_fingerprint = {}
    for pubchem_id, smiles in tqdm(node_to_smiles.items(), desc='Getting fingerprints'):
        mol_from_smile = Chem.MolFromSmiles(smiles)
        if mol_from_smile is None:
            continue
        node_to_fingerprint[pubchem_id] = MACCSkeys.GenMACCSKeys(mol_from_smile)
    return node_to_fingerprint
