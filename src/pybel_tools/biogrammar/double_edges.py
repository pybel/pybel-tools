# -*- coding: utf-8 -*-

"""This module investigates the properties of paths of length two (A - B - C)."""

from typing import Dict, Iterable, List, Mapping, Tuple

from pybel import BELGraph
from pybel.constants import (
    ACTIVITY, CAUSAL_RELATIONS, DEGRADATION, HAS_COMPONENT, LINE, MODIFIER, OBJECT, RELATION, SUBJECT, TRANSLOCATION,
)
from pybel.dsl import BaseEntity, ComplexAbundance, NamedComplexAbundance
from pybel.struct.filters import edge_predicate, has_protein_modification, part_has_modifier
from pybel.typing import EdgeData


def self_edge_filter(_: BELGraph, source: BaseEntity, target: BaseEntity, __: str) -> bool:
    """Check if the source and target nodes are the same."""
    return source == target


def has_protein_modification_increases_activity(
        graph: BELGraph,
        source: BaseEntity,
        target: BaseEntity,
        key: str,
) -> bool:
    """Check if pmod of source causes activity of target."""
    edge_data = graph[source][target][key]
    return has_protein_modification(graph, source) and part_has_modifier(edge_data, OBJECT, ACTIVITY)


@edge_predicate
def has_degradation_increases_activity(edge_data: EdgeData) -> bool:
    """Check if the degradation of source causes activity of target."""
    return part_has_modifier(edge_data, SUBJECT, DEGRADATION) and part_has_modifier(edge_data, OBJECT, ACTIVITY)


@edge_predicate
def has_translocation_increases_activity(edge_data: EdgeData) -> bool:
    """Check if the translocation of source causes activity of target."""
    return part_has_modifier(edge_data, SUBJECT, TRANSLOCATION) and part_has_modifier(edge_data, OBJECT, ACTIVITY)


def complex_has_member(graph: BELGraph, complex_node: ComplexAbundance, member_node: BaseEntity) -> bool:
    """Check if the given complex contains the member."""
    return any(  # TODO can't you look in the members of the complex object (if it's enumerated)
        v == member_node
        for _, v, data in graph.out_edges(complex_node, data=True)
        if data[RELATION] == HAS_COMPONENT
    )


def complex_increases_activity(graph: BELGraph, u: BaseEntity, v: BaseEntity, key: str) -> bool:
    """Return if the formation of a complex with u increases the activity of v."""
    return (
        isinstance(u, (ComplexAbundance, NamedComplexAbundance)) and
        complex_has_member(graph, u, v) and
        part_has_modifier(graph[u][v][key], OBJECT, ACTIVITY)
    )


def has_same_subject_object(graph: BELGraph, u: BaseEntity, v: BaseEntity, key: str) -> bool:
    """Check if the edge has the same subject and object modifiers."""
    data = graph[u][v][key]
    return u == v and data.get(SUBJECT) == data.get(OBJECT)


def get_related_causal_out_edges(
        graph: BELGraph,
        node: BaseEntity,
        edge_data: EdgeData,
) -> Iterable[Tuple[BaseEntity, Dict, str]]:
    """Get downstream nodes whose activity is modulated by the given node."""
    for _, v, data in graph.out_edges(node, data=True):
        if data[RELATION] in CAUSAL_RELATIONS and data.get(SUBJECT) == edge_data:
            yield v, data, graph.edge_to_bel(node, v, data)


def get_related_ish_causal_out_edges(
        graph: BELGraph,
        node: BaseEntity,
        modifier: str,
) -> Iterable[Tuple[BaseEntity, Dict, str]]:
    """Get targets and relations to the given node that are modified."""
    for _, v, data in graph.out_edges(node, data=True):
        if data[RELATION] in CAUSAL_RELATIONS and data.get(SUBJECT, {}).get(MODIFIER) == modifier:
            yield v, data, graph.edge_to_bel(node, v, data)


def summarize_completeness(graph: BELGraph) -> Mapping[str, List]:
    """Summarize the completeness of a BEL graph using several categories.

    1. Protein modifications that result in a change in activity
    2. Degradations that result in a change in activity
    3. Transloations that result in a change in activity
    4. Complex formation that results in a change in activity
    """
    pmod_activity = []
    deg_activity = []
    tloc_activity = []
    comp_activity = []
    same_sub_obj = []
    same_undefined = []

    for u, v, key, data in graph.edges(keys=True, data=True):
        if u != v:
            continue

        bel = graph.edge_to_bel(u, v, data)

        line = data.get(LINE)

        if line is None:
            continue  # this was inferred, so need to investigate another way

        entry = {
            'bel': bel,
            'line': line,
        }

        if has_protein_modification_increases_activity(graph, u, v, key):
            e = entry.copy()
            e['chain'] = [
                {
                    'bel': v_w_bel,
                }
                for w, v_w_data, v_w_bel in get_related_causal_out_edges(graph, v, data)
            ]
            e['related'] = [
                {
                    'bel': v_w_bel,
                }
                for w, v_w_data, v_w_bel in get_related_ish_causal_out_edges(graph, v, data)
            ]
            pmod_activity.append(e)

        elif has_degradation_increases_activity(data):
            e = entry.copy()
            e['chain'] = [
                {
                    'bel': v_w_bel,
                }
                for w, v_w_data, v_w_bel in get_related_causal_out_edges(graph, v, data)
            ]
            e['related'] = [
                {
                    'bel': v_w_bel,
                }
                for w, v_w_data, v_w_bel in get_related_ish_causal_out_edges(graph, v, data)
            ]
            deg_activity.append(e)

        elif has_translocation_increases_activity(data):
            e = entry.copy()
            e['chain'] = [
                {
                    'bel': v_w_bel,
                }
                for w, v_w_data, v_w_bel in get_related_causal_out_edges(graph, v, data)
            ]
            e['related'] = [
                {
                    'bel': v_w_bel,
                }
                for w, v_w_data, v_w_bel in get_related_ish_causal_out_edges(graph, v, data)
            ]
            tloc_activity.append(e)

        elif complex_increases_activity(graph, u, v, key):
            e = entry.copy()
            e['chain'] = [
                {
                    'bel': v_w_bel,
                }
                for w, v_w_data, v_w_bel in get_related_causal_out_edges(graph, v, data)
            ]
            e['related'] = [
                {
                    'bel': v_w_bel,
                }
                for w, v_w_data, v_w_bel in get_related_ish_causal_out_edges(graph, v, data)
            ]
            comp_activity.append(e)

        elif has_same_subject_object(graph, u, v, key):
            same_sub_obj.append(entry.copy())

        else:
            same_undefined.append(entry.copy())

    return dict(
        pmod_activity=pmod_activity,
        deg_activity=deg_activity,
        tloc_activity=tloc_activity,
        comp_activity=comp_activity,
        same_sub_obj=same_sub_obj,
        same_undefined=same_undefined,
    )


def find_activations(graph: BELGraph):
    """Find edges that are A - A, meaning that some conditions in the edge best describe the interaction."""
    for u, v, key, data in graph.edges(keys=True, data=True):
        if u != v:
            continue

        bel = graph.edge_to_bel(u, v, data)

        line = data.get(LINE)

        if line is None:
            continue  # this was inferred, so need to investigate another way

        elif has_protein_modification_increases_activity(graph, u, v, key):
            print(line, '- pmod changes -', bel)
            find_related(graph, v, data)

        elif has_degradation_increases_activity(data):
            print(line, '- degradation changes -', bel)
            find_related(graph, v, data)

        elif has_translocation_increases_activity(data):
            print(line, '- translocation changes -', bel)
            find_related(graph, v, data)

        elif complex_increases_activity(graph, u, v, key):
            print(line, '- complex changes - ', bel)
            find_related(graph, v, data)

        elif has_same_subject_object(graph, u, v, key):
            print(line, '- same sub/obj -', bel)

        else:
            print(line, '- *** - ', bel)


def find_related(graph: BELGraph, v: BaseEntity, data):
    """Find related edges."""
    for bel in sorted({bel for _, _, bel in get_related_causal_out_edges(graph, v, data)}):
        print('* chained - ', bel)

    for bel in sorted({bel for _, _, bel in get_related_ish_causal_out_edges(graph, v, data.get(MODIFIER))}):
        print('^ related - ', bel)
