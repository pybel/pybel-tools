# -*- coding: utf-8 -*-

"""Generate HTML summaries of BEL graphs."""

from __future__ import annotations

import itertools as itt
import logging
import os
from dataclasses import dataclass
from typing import Counter, List, Mapping, Optional, TextIO, Tuple

from dataclasses_json import dataclass_json
from jinja2 import Environment, FileSystemLoader

from pybel import BELGraph
from pybel.constants import IDENTIFIER, NAME, RELATION
from pybel.dsl import BaseEntity
from pybel.struct.filters import filter_edges
from pybel.struct.summary import (
    count_functions, count_namespaces, count_relations, count_variants, get_syntax_errors, get_top_hubs,
    get_top_pathologies, get_unused_namespaces,
)
from ...analysis.stability import (
    get_chaotic_pairs, get_chaotic_triplets, get_contradiction_summary, get_dampened_pairs, get_dampened_triplets,
    get_decrease_mismatch_triplets, get_increase_mismatch_triplets, get_jens_unstable,
    get_mutually_unstable_correlation_triples, get_regulatory_pairs, get_separate_unstable_correlation_triples,
)
from ...filters import has_pathology_causal
from ...summary import (
    count_authors, count_confidences, count_error_types, get_citation_years, get_modifications_count,
    get_most_common_errors, get_naked_names, get_namespaces_with_incorrect_names, get_undefined_annotations,
    get_undefined_namespaces, get_unused_annotations, get_unused_list_annotation_values,
)
from ...typing import SetOfNodePairs, SetOfNodeTriples
from ...utils import prepare_c3, prepare_c3_time_series

__all__ = [
    'to_html',
    'to_html_file',
    'to_html_path',
    'get_network_summary_dict',
    'prepare_c3',
]

log = logging.getLogger(__name__)

HERE = os.path.dirname(os.path.abspath(__file__))
environment = Environment(autoescape=False, loader=FileSystemLoader(HERE), trim_blocks=False)


def to_html_path(graph: BELGraph, path: str) -> None:
    """Write the graph as HTML to a file at the given path."""
    with open(path, 'w') as file:
        to_html_file(graph, file)


def to_html_file(graph: BELGraph, file: Optional[TextIO] = None) -> None:
    """Write the graph as HTML to a file."""
    html = to_html(graph)
    print(html, file=file)


def to_html(graph: BELGraph) -> str:
    """Render the graph as an HTML string.

    Common usage may involve writing to a file like:

    >>> from pybel.examples import sialic_acid_graph
    >>> with open('html_output.html', 'w') as file:
    ...     print(to_html(sialic_acid_graph), file=file)
    """
    context = get_network_summary_dict(graph)
    summary_dict = graph.summary_dict()

    citation_years = context['citation_years']
    function_count = context['function_count']
    relation_count = context['relation_count']
    error_count = context['error_count']
    transformations_count = context['modifications_count']
    hub_data = context['hub_data']
    disease_data = context['disease_data']
    authors_count = context['authors_count']
    variants_count = context['variants_count']
    namespaces_count = context['namespaces_count']
    confidence_count = context['confidence_count']
    confidence_data = [
        (label, confidence_count.get(label, 0))
        for label in ('None', 'Very Low', 'Low', 'Medium', 'High', 'Very High')
    ]

    template = environment.get_template('index.html')
    return template.render(
        graph=graph,
        # Node Charts
        chart_1_data=prepare_c3(function_count, 'Node Count'),
        chart_6_data=prepare_c3(namespaces_count, 'Node Count'),
        chart_5_data=prepare_c3(variants_count, 'Node Count'),
        number_variants=sum(variants_count.values()),
        number_namespaces=len(namespaces_count),
        # Edge Charts
        chart_2_data=prepare_c3(relation_count, 'Edge Count'),
        chart_4_data=prepare_c3(transformations_count, 'Edge Count') if transformations_count else None,
        number_transformations=sum(transformations_count.values()),
        # Error Charts
        chart_3_data=prepare_c3(error_count, 'Error Type') if error_count else None,
        # Topology Charts
        chart_7_data=prepare_c3(hub_data, 'Degree'),
        chart_9_data=prepare_c3(disease_data, 'Degree') if disease_data else None,
        # Bibliometrics Charts
        chart_authors_count=prepare_c3(authors_count.most_common(15), 'Edges Contributed'),
        chart_10_data=prepare_c3_time_series(citation_years, 'Number of Articles') if citation_years else None,
        chart_confidence_count=prepare_c3(confidence_data, 'Edge Count'),
        summary_dict=summary_dict,
        # Everything else :)
        **context
    )


@dataclass_json
@dataclass
class BELGraphSummary:
    """A container for a summary of a BEL graph."""

    # Attribute counters
    function_count: Counter[str]
    modifications_count: Counter[str]
    relation_count: Counter[str]
    authors_count: Counter[str]
    variants_count: Counter[str]
    namespaces_count: Counter[str]

    # Node counters
    hubs: Counter[BaseEntity]

    # Node pairs
    regulatory_pairs: SetOfNodePairs
    chaotic_pairs: SetOfNodePairs
    dampened_pairs: SetOfNodePairs
    contradictory_pairs: SetOfNodePairs

    # Node triplets
    separate_unstable_correlation_triples: SetOfNodeTriples
    mutually_unstable_correlation_triples: SetOfNodeTriples
    jens_unstable: SetOfNodeTriples
    increase_mismatch_triplets: SetOfNodeTriples
    decrease_mismatch_triplets: SetOfNodeTriples

    @staticmethod
    def from_graph(graph: BELGraph) -> BELGraphSummary:
        """Create a summary of the graph."""
        return BELGraphSummary(
            function_count=count_functions(graph),
            modifications_count=get_modifications_count(graph),
            relation_count=count_relations(graph),
            authors_count=count_authors(graph),
            variants_count=count_variants(graph),
            namespaces_count=count_namespaces(graph),
            regulatory_pairs=get_regulatory_pairs(graph),
            chaotic_pairs=get_chaotic_pairs(graph),
            dampened_pairs=get_dampened_pairs(graph),
            contradictory_pairs=get_contradiction_summary(graph),
            separate_unstable_correlation_triples=get_separate_unstable_correlation_triples(graph),
            mutually_unstable_correlation_triples=get_mutually_unstable_correlation_triples(graph),
            jens_unstable=get_jens_unstable(graph),
            increase_mismatch_triplets=get_increase_mismatch_triplets(graph),
            decrease_mismatch_triplets=get_decrease_mismatch_triplets(graph),
        )


def get_network_summary_dict(graph: BELGraph) -> Mapping:
    """Create a summary dictionary."""
    summary = BELGraphSummary.from_graph(graph)

    contradictory_triplets = list(itt.chain(
        (get_triplet_tuple(a, b, c) + ('Separate',) for a, b, c in summary.separate_unstable_correlation_triples),
        (get_triplet_tuple(a, b, c) + ('Mutual',) for a, b, c in summary.mutually_unstable_correlation_triples),
        (get_triplet_tuple(a, b, c) + ('Jens',) for a, b, c in summary.jens_unstable),
        (get_triplet_tuple(a, b, c) + ('Increase Mismatch',) for a, b, c in summary.increase_mismatch_triplets),
        (get_triplet_tuple(a, b, c) + ('Decrease Mismatch',) for a, b, c in summary.decrease_mismatch_triplets),
    )),

    regulatory_pairs = [
        get_pair_tuple(u, v)
        for u, v in summary.regulatory_pairs
    ]

    unstable_pairs = list(itt.chain(
        (get_pair_tuple(u, v) + ('Chaotic',) for u, v in summary.chaotic_pairs),
        (get_pair_tuple(u, v) + ('Dampened',) for u, v in get_dampened_pairs(graph)),
    ))

    contradictory_pairs = [
        get_pair_tuple(u, v) + (relation,)
        for u, v, relation in summary.contradictory_pairs
    ]

    rv = dict(
        undefined_namespaces=get_undefined_namespaces(graph),
        undefined_annotations=get_undefined_annotations(graph),
        namespaces_with_incorrect_names=get_namespaces_with_incorrect_names(graph),
        unused_namespaces=get_unused_namespaces(graph),
        unused_annotations=get_unused_annotations(graph),
        unused_list_annotation_values=get_unused_list_annotation_values(graph),
        naked_names=get_naked_names(graph),
        error_count=count_error_types(graph),
        # Errors
        error_groups=get_most_common_errors(graph),
        syntax_errors=get_syntax_errors(graph),
        # Bibliometrics
        citation_years=get_citation_years(graph),
        confidence_count=count_confidences(graph),

        causal_pathologies=sorted({
            get_pair_tuple(u, v) + (graph[u][v][k][RELATION],)
            for u, v, k in filter_edges(graph, has_pathology_causal)
        }),
    )

    rv['hub_data'] = {
        (
            node.name or node.identifier
            if NAME in node or IDENTIFIER in node else
            str(node)
        ): degree
        for node, degree in get_top_hubs(graph, n=15)
    }
    rv['disease_data'] = {
        (
            node.name or node.identifier
            if NAME in node or IDENTIFIER in node else
            str(node)
        ): count
        for node, count in get_top_pathologies(graph, n=15)
    }

    rv.update(summary.to_dict())

    return rv


BELPairTuple = Tuple[str, str, str, str]
BELTripleTuple = Tuple[str, str, str, str, str, str]


def get_pair_tuple(a: BaseEntity, b: BaseEntity) -> BELPairTuple:
    """Get the pair as a tuple of BEL/hashes."""
    return a.as_bel(), a.sha512, b.as_bel(), b.sha512


def get_triplet_tuple(a: BaseEntity, b: BaseEntity, c: BaseEntity) -> BELTripleTuple:
    """Get the triple as a tuple of BEL/hashes."""
    return a.as_bel(), a.sha512, b.as_bel(), b.sha512, c.as_bel(), c.sha512
