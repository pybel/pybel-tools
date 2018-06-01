# -*- coding: utf-8 -*-

"""This module investigates the properties of paths of length two (A - B - C)"""

from pybel.constants import (
    ACTIVITY, CAUSAL_RELATIONS, COMPLEX, DEGRADATION, FUNCTION, HAS_COMPONENT, LINE, MODIFIER, OBJECT, RELATION,
    SUBJECT, TRANSLOCATION,
)
from pybel.struct.filters import edge_predicate, has_protein_modification, part_has_modifier


def self_edge_filter(graph, source, target, key):
    """Checks if the source and target nodes are the same

    :param pybel.BELGraph graph:
    :param tuple source:
    :param tuple target:
    :param int key:
    :rtype: bool
    """
    return source == target


def has_protein_modification_increases_activity(graph, source, target, key):
    """Check if pmod of source causes activity of target

    :param graph:
    :param source:
    :param target:
    :param key:
    :rtype:
    """
    edge_data = graph.edge[source][target][key]

    return has_protein_modification(graph, source) and part_has_modifier(edge_data, OBJECT, ACTIVITY)


@edge_predicate
def has_degradation_increases_activity(data):
    """Check if the degradation of source causes activity of target

    :param dict data:
    """
    return part_has_modifier(data, SUBJECT, DEGRADATION) and part_has_modifier(data, OBJECT, ACTIVITY)


@edge_predicate
def has_translocation_increases_activity(data):
    """Check if the tranclocation of source causes activity of target

    :param dict data:
    """
    return part_has_modifier(data, SUBJECT, TRANSLOCATION) and part_has_modifier(data, OBJECT, ACTIVITY)


def complex_has_member(graph, complex_node, member_node):
    """Does the given complex contain the member?

    :param pybel.BELGraph graph:
    :param tuple complex_node:
    :param tuple member_node:
    :rtype: bool
    """
    for _, v, data in graph.out_edges(complex_node, data=True):
        if data[RELATION] != HAS_COMPONENT:
            continue
        if v == member_node:
            return True


def complex_increases_activity(graph, u, v, key):
    """If the complexing of v increases its activity"""
    if graph.node[u][FUNCTION] != COMPLEX:
        return False

    if not complex_has_member(graph, u, v):
        return False

    return part_has_modifier(graph.edge[u][v][key], OBJECT, ACTIVITY)


def has_same_subject_object(graph, u, v, key):
    data = graph.edge[u][v][key]
    return u == v and data.get(SUBJECT) == data.get(OBJECT)


def get_related_causal_out_edges(graph, u, edge_data):
    for _, v, data in graph.out_edges(u, data=True):
        if data[RELATION] in CAUSAL_RELATIONS and data.get(SUBJECT) == edge_data:
            yield v, data, graph.edge_to_bel(u, v, data)


def get_related_ish_causal_out_edges(graph, u, modifier):
    for _, v, data in graph.out_edges(u, data=True):
        if data[RELATION] in CAUSAL_RELATIONS and data.get(SUBJECT, {}).get(MODIFIER) == modifier:
            yield v, data, graph.edge_to_bel(u, v, data)


def find_related(graph, v, data):
    for bel in sorted({bel for _, _, bel in get_related_causal_out_edges(graph, v, data)}):
        print('* chained - ', bel)

    for bel in sorted({bel for _, _, bel in get_related_ish_causal_out_edges(graph, v, data.get(MODIFIER))}):
        print('^ related - ', bel)


def summarize_competeness(graph):
    """

    :param pybel.BELGraph graph:
    :rtype: list[dict]
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
        same_undefined=same_undefined
    )


def find_activations(graph):
    """Find edges that are A - A, meaning that some conditions in the edge best describe the interaction

    :param pybel.BELGraph graph: A BEL graph
    """
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


if __name__ == '__main__':
    import pybel

    # g = pybel.from_pickle('/Users/cthoyt/dev/bms/aetionomy/alzheimers/alzheimers.gpickle')
    g = pybel.from_pickle('/Users/cthoyt/dev/bms/selventa/large_corpus.gpickle')
    find_activations(g)
