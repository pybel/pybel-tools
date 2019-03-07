# -*- coding: utf-8 -*-

"""Utilities for handling nodes."""

import itertools as itt
import logging
from itertools import chain

from networkx import relabel_nodes
from pybel import BELGraph
from pybel.constants import ANNOTATIONS, CITATION, EVIDENCE, RELATION
from pybel.dsl import ListAbundance, Reaction

__all__ = [
    'flatten_list_abundance',
    'list_abundance_cartesian_expansion',
    'reaction_cartesian_expansion',
]

logger = logging.getLogger(__name__)


def flatten_list_abundance(node: ListAbundance) -> ListAbundance:
    """Flattens the complex or composite abundance."""
    return node.__class__(list(chain.from_iterable(
        (
            flatten_list_abundance(member).members
            if isinstance(member, ListAbundance) else
            [member]
        )
        for member in node.members
    )))


def list_abundance_expansion(graph: BELGraph) -> None:
    """Flatten list abundances"""
    mapping = {
        node: flatten_list_abundance(node)
        for node in graph
        if isinstance(node, ListAbundance)
    }
    relabel_nodes(graph, mapping, copy=False)


def list_abundance_cartesian_expansion(graph: BELGraph) -> None:
    """Expand all list abundances to simple subject-predicate-object networks."""
    for u, v, k, d in list(graph.edges(keys=True, data=True)):
        if CITATION not in d:
            continue

        if isinstance(u, ListAbundance) and isinstance(v, ListAbundance):
            for u_member, v_member in itt.product(u.members, v.members):
                graph.add_qualified_edge(
                    u_member, v_member,
                    relation=d[RELATION],
                    citation=d.get(CITATION),
                    evidence=d.get(EVIDENCE),
                    annotations=d.get(ANNOTATIONS),
                )

        elif isinstance(u, ListAbundance):
            for member in u.members:
                graph.add_qualified_edge(
                    member, v,
                    relation=d[RELATION],
                    citation=d.get(CITATION),
                    evidence=d.get(EVIDENCE),
                    annotations=d.get(ANNOTATIONS),
                )

        elif isinstance(v, ListAbundance):
            for member in v.members:
                graph.add_qualified_edge(
                    u, member,
                    relation=d[RELATION],
                    citation=d.get(CITATION),
                    evidence=d.get(EVIDENCE),
                    annotations=d.get(ANNOTATIONS),
                )

    list_nodes = {
        node
        for node in graph.nodes()
        if isinstance(node, ListAbundance)
    }
    graph.remove_nodes_from(list_nodes)


def reaction_cartesian_expansion(graph: BELGraph) -> None:
    """Expand all reactions to simple subject-predicate-object networks."""
    for u, v, k, d in list(graph.edges(keys=True, data=True)):
        if CITATION not in d:
            continue

        if isinstance(u, Reaction) and isinstance(v, Reaction):
            for reactant, product in chain(itt.product(u.reactants, u.products), itt.product(v.reactants, v.products)):
                graph.add_increases(
                    reactant, product,
                    citation=d.get(CITATION),
                    evidence=d.get(EVIDENCE),
                    annotations=d.get(ANNOTATIONS),
                )
            for product, reactant in itt.product(u.products, u.reactants):
                graph.add_qualified_edge(
                    product, reactant,
                    relation=d[RELATION],
                    citation=d.get(CITATION),
                    evidence=d.get(EVIDENCE),
                    annotations=d.get(ANNOTATIONS),
                )

        elif isinstance(u, Reaction):
            for product in u.products:
                graph.add_increases(
                    product, v,
                    citation=d.get(CITATION),
                    evidence=d.get(EVIDENCE),
                    annotations=d.get(ANNOTATIONS),
                )
                for reactant in u.reactants:
                    graph.add_increases(
                        reactant, product,
                        citation=d.get(CITATION),
                        evidence=d.get(EVIDENCE),
                        annotations=d.get(ANNOTATIONS),
                    )

        elif isinstance(v, Reaction):
            for reactant in v.reactants:
                graph.add_increases(
                    u, reactant,
                    citation=d.get(CITATION),
                    evidence=d.get(EVIDENCE),
                    annotations=d.get(ANNOTATIONS),
                )
                for product in v.products:
                    graph.add_increases(
                        reactant, product,
                        citation=d.get(CITATION),
                        evidence=d.get(EVIDENCE),
                        annotations=d.get(ANNOTATIONS),
                    )

    reaction_nodes = {
        node
        for node in graph.nodes()
        if isinstance(node, Reaction)
    }
    graph.remove_nodes_from(reaction_nodes)

def remove_complex_nodes(graph:BELGraph) -> None:
    """Remove complex nodes."""
    list_nodes = {
        node
        for node in graph.nodes()
        if isinstance(node, ListAbundance)
    }
    graph.remove_nodes_from(list_nodes)

    reaction_nodes = {
        node
        for node in graph.nodes()
        if isinstance(node, Reaction)
    }
    graph.remove_nodes_from(reaction_nodes)
