# -*- coding: utf-8 -*-

"""Utilities for handling nodes."""

import itertools as itt
import logging
from itertools import chain

from networkx import relabel_nodes
from pybel import BELGraph
from pybel.constants import ANNOTATIONS, CITATION, EVIDENCE, RELATION, INCREASES, HAS_REACTANT
from pybel.dsl import BaseEntity, ListAbundance, Reaction

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
    """Flatten list abundances."""
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


def _reaction_cartesion_expansion_unqualified_helper(graph: BELGraph, u: BaseEntity, v: BaseEntity, d: dict) -> None:
    """Helper to deal with cartension expansion in unqualified edges."""
    if isinstance(u, Reaction) and isinstance(v, Reaction):

        enzymes = _get_enzymes_in_reaction(u) + _get_enzymes_in_reaction(v)

        for reactant, product in chain(itt.product(u.reactants, u.products),
                                       itt.product(v.reactants, v.products)):

            if reactant in enzymes or product in enzymes:
                continue

            graph.add_unqualified_edge(
                reactant, product, INCREASES
            )
        for product, reactant in itt.product(u.products, u.reactants):

            if reactant in enzymes or product in enzymes:
                continue

            graph.add_unqualified_edge(
                product, reactant, d[RELATION],
            )

    elif isinstance(u, Reaction):

        enzymes = _get_enzymes_in_reaction(u)

        for product in u.products:

            # Skip create increases edges between enzymes
            if product in enzymes:
                continue

            # Only add edge between v and reaction if the node is not part of the reaction
            # In practice skips hasReactant, hasProduct edges
            if v not in u.products and v not in u.reactants:
                graph.add_unqualified_edge(
                    product, v, INCREASES
                )
            for reactant in u.reactants:
                graph.add_unqualified_edge(
                    reactant, product, INCREASES
                )

    elif isinstance(v, Reaction):

        enzymes = _get_enzymes_in_reaction(v)

        for reactant in v.reactants:

            # Skip create increases edges between enzymes
            if reactant in enzymes:
                continue

            # Only add edge between v and reaction if the node is not part of the reaction
            # In practice skips hasReactant, hasProduct edges
            if u not in v.products and u not in v.reactants:
                graph.add_unqualified_edge(
                    u, reactant, INCREASES
                )
            for product in v.products:
                graph.add_unqualified_edge(
                    reactant, product, INCREASES
                )


def _get_enzymes_in_reaction(reaction: Reaction):
    """Return nodes that are both in reactants and reactions in a reaction."""
    return [
        reactant
        for reactant in reaction.reactants
        if reactant in reaction.products
    ]


def reaction_cartesian_expansion(graph: BELGraph, accept_unqualified_edges=True) -> None:
    """Expand all reactions to simple subject-predicate-object networks."""
    for u, v, k, d in list(graph.edges(keys=True, data=True)):

        # Deal with unqualified edges
        if CITATION not in d and accept_unqualified_edges:
            _reaction_cartesion_expansion_unqualified_helper(graph, u, v, d)
            continue

        if isinstance(u, Reaction) and isinstance(v, Reaction):

            enzymes = _get_enzymes_in_reaction(u) + _get_enzymes_in_reaction(v)

            for reactant, product in chain(itt.product(u.reactants, u.products), itt.product(v.reactants, v.products)):

                if reactant in enzymes or product in enzymes:
                    continue

                graph.add_increases(
                    reactant, product,
                    citation=d.get(CITATION),
                    evidence=d.get(EVIDENCE),
                    annotations=d.get(ANNOTATIONS),
                )
            for product, reactant in itt.product(u.products, u.reactants):

                if reactant in enzymes or product in enzymes:
                    continue

                graph.add_qualified_edge(
                    product, reactant,
                    relation=d[RELATION],
                    citation=d.get(CITATION),
                    evidence=d.get(EVIDENCE),
                    annotations=d.get(ANNOTATIONS),
                )

        elif isinstance(u, Reaction):

            enzymes = _get_enzymes_in_reaction(u)

            for product in u.products:

                # Skip create increases edges between enzymes
                if product in enzymes:
                    continue

                # Only add edge between v and reaction if the node is not part of the reaction
                # In practice skips hasReactant, hasProduct edges
                if v not in u.products and v not in u.reactants:
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
                enzymes = _get_enzymes_in_reaction(v)

                # Skip create increases edges between enzymes
                if reactant in enzymes:
                    continue

                # Only add edge between v and reaction if the node is not part of the reaction
                # In practice skips hasReactant, hasProduct edges
                if u not in v.products and u not in v.reactants:
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


def remove_complex_nodes(graph: BELGraph) -> None:
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
