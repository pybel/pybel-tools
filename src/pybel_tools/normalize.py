# -*- coding: utf-8 -*-

"""Utilities to help with normalizing entities in BEL graphs."""

import logging
from collections import Counter
from itertools import chain
from typing import Iterable, List, Mapping, Optional, TextIO, Tuple

from pybel import BELGraph
from pybel.constants import IDENTIFIER, NAME, NAMESPACE
from pybel.dsl import (
    BaseEntity, CentralDogma, FusionBase, GeneModification, ListAbundance,
    ProteinModification, Reaction,
)

from pybel_tools.utils import group_as_lists

__all__ = [
    'normalize',
    'get_unnormalized',
    'summarize_unnormalized',
]

logger = logging.getLogger(__name__)


def normalize(graph: BELGraph, use_tqdm: bool = True) -> None:
    """Normalize all of the entities in the graph."""
    logger.info('normalizing HGNC')
    import bio2bel_hgnc
    hgnc_manager = bio2bel_hgnc.Manager()
    hgnc_manager.normalize_genes(graph, use_tqdm=use_tqdm)

    # logger.info('normalizing HGNC Gene Families')
    # gfam_manager = bio2bel_hgnc.FamilyManager()
    # gfam_manager.normalize_families(graph)

    logger.info('normalizing MeSH')
    import bio2bel_mesh
    mesh_manager = bio2bel_mesh.Manager()
    if mesh_manager.is_populated():
        mesh_manager.normalize_terms(graph, use_tqdm=use_tqdm)
    else:
        logger.warning('MeSH has not been populated')

    logger.info('normalizing FamPlex')
    import bio2bel_famplex
    famplex_manager = bio2bel_famplex.Manager()
    famplex_manager.normalize_terms(graph, use_tqdm=use_tqdm)

    logger.info('normalizing GO')
    import bio2bel_go
    go_manager = bio2bel_go.Manager()
    go_manager.normalize_terms(graph, use_tqdm=use_tqdm)

    logger.info('normalizing CONSO')
    import conso.manager
    conso_manager = conso.manager.Manager()
    conso_manager.normalize_terms(graph, use_tqdm=use_tqdm)

    logger.info('normalizing ChEBI')
    import bio2bel_chebi
    chebi_manager = bio2bel_chebi.Manager()
    if chebi_manager.is_populated():
        chebi_manager.normalize_chemicals(graph, use_tqdm=use_tqdm)
    else:
        logger.warning('ChEBI has not been populated')

    logger.info('normalizing MGI')
    import bio2bel_mgi
    mgi_manager = bio2bel_mgi.Manager()
    if mgi_manager.is_populated():
        mgi_manager.normalize_mouse_genes(graph, use_tqdm=use_tqdm)
    else:
        logger.warning('MGI has not been populated')

    logger.info('normalizing RGD')
    import bio2bel_rgd
    rgd_manager = bio2bel_rgd.Manager()
    if rgd_manager.is_populated():
        rgd_manager.normalize_rat_genes(graph, use_tqdm=use_tqdm)
    else:
        logger.warning('RGD has not been populated')

    # logger.info('normalizing InterPro')

    # logger.info('normalizing PFAM')

    logger.info('normalizing Entrez Gene')
    import bio2bel_entrez
    entrez_manager = bio2bel_entrez.Manager()
    entrez_manager.normalize_genes(graph, use_tqdm=use_tqdm)

    logger.info('normalizing UniProt')
    import bio2bel_uniprot
    uniprot_manager = bio2bel_uniprot.Manager()
    uniprot_manager.normalize_terms(graph, use_tqdm=use_tqdm)

    logger.info('normalizing DrugBank')
    import bio2bel_drugbank
    drugbank_manager = bio2bel_drugbank.Manager()
    if drugbank_manager.is_populated():
        drugbank_manager.normalize_drugs(graph, use_tqdm=use_tqdm)
    else:
        logger.warning('DrugBank has not been populated')

    logger.info('normalizing miRBase')
    import bio2bel_mirbase
    mirbase_manager = bio2bel_mirbase.Manager()
    if mirbase_manager.is_populated():
        mirbase_manager.normalize_terms(graph, use_tqdm=use_tqdm)
    else:
        logger.warning('miRBase has not been populated')

    # TODO deal with OBO-based bio2bel
    # logger.info('normalizing CL')
    # logger.info('normalizing HP')


def summarize_unnormalized(graph: BELGraph, file: Optional[TextIO] = None) -> None:
    for namespace, names in get_unnormalized(graph).items():
        name_counter = Counter(names)
        print(
            namespace,
            len(names),
            *[x for x, _ in name_counter.most_common(3)],
            file=file,
        )


def get_unnormalized(graph: BELGraph) -> Mapping[str, List[str]]:
    """Get a mapping of namespaces to their unnormalized names."""
    return group_as_lists(iter_unnormalize_entity_namespaces(graph))


def iter_unnormalize_entity_namespaces(graph: BELGraph) -> Iterable[Tuple[str, str]]:
    """Get the namespaces that haven't been normalized."""
    for node in graph:
        yield from _iter_unnormalized_node(node)


def _iter_unnormalized_node(node: BaseEntity) -> Iterable[Tuple[str, str]]:
    namespace, name, identifier = node.get(NAMESPACE), node.get(NAME), node.get(IDENTIFIER)

    if not namespace:
        pass

    elif name and identifier:
        pass

    elif name and not identifier:
        yield namespace, name

    elif not name and identifier:
        yield namespace, identifier

    elif isinstance(node, FusionBase):
        yield from _iter_unnormalized_node(node.partner_5p)
        yield from _iter_unnormalized_node(node.partner_3p)

    elif isinstance(node, CentralDogma) and node.variants:
        for variant in node.variants:
            if isinstance(variant, (GeneModification, ProteinModification)):
                namespace = variant.entity.get(NAMESPACE)
                name = variant.entity.get(NAME)
                identifier = variant.entity.get(IDENTIFIER)
                if namespace and not identifier:
                    yield namespace, name

    elif isinstance(node, ListAbundance):
        for member in node.members:
            yield from _iter_unnormalized_node(member)

    elif isinstance(node, Reaction):
        for member in chain(node.reactants, node.products):
            yield from _iter_unnormalized_node(member)

    else:
        logger.warning('Unhandled node: %r', node)
