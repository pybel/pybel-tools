# -*- coding: utf-8 -*-

"""A BEL Graph summary class."""

from __future__ import annotations

import collections
from dataclasses import dataclass
from typing import Counter, List, Mapping, Set, Tuple

from dataclasses_json import dataclass_json

from pybel import BELGraph, BaseAbundance, BaseEntity
from pybel.constants import IDENTIFIER, NAME
from pybel.struct.graph import WarningTuple
from pybel.struct.summary import (
    count_error_types, count_functions, count_namespaces, count_relations, count_variants, get_naked_names,
    get_syntax_errors, get_top_hubs, get_top_pathologies, get_unused_annotations, get_unused_list_annotation_values,
    get_unused_namespaces,
)
from .error_summary import (
    get_most_common_errors, get_namespaces_with_incorrect_names, get_undefined_annotations,
    get_undefined_namespaces,
)
from .node_properties import count_modifications
from .provenance import count_authors, count_confidences, get_citation_years
from .stability import (
    get_chaotic_pairs, get_contradiction_summary, get_dampened_pairs, get_decrease_mismatch_triplets,
    get_increase_mismatch_triplets, get_jens_unstable, get_mutually_unstable_correlation_triples, get_regulatory_pairs,
    get_separate_unstable_correlation_triples,
)
from ..typing import SetOfNodePairs, SetOfNodeTriples

__all__ = [
    'BELGraphSummary',
]


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

    # Errors
    undefined_namespaces: Set[str]
    undefined_annotations: Set[str]
    namespaces_with_incorrect_names: Set[str]
    unused_namespaces: Set[str]
    unused_annotations: Set[str]
    unused_list_annotation_values: Mapping[str, Set[str]]
    naked_names: Set[str]
    error_count: Counter[str]
    error_groups: List[Tuple[str, str]]
    syntax_errors: List[WarningTuple]

    # Node counters
    hub_data: Counter[BaseEntity]
    disease_data: Counter[BaseEntity]

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

    # Bibliometrics
    citation_years: List[Tuple[int, int]]
    confidence_count: Counter[str]

    @staticmethod
    def from_graph(graph: BELGraph) -> BELGraphSummary:
        """Create a summary of the graph."""
        return BELGraphSummary(
            # Attribute counters
            function_count=count_functions(graph),
            modifications_count=count_modifications(graph),
            relation_count=count_relations(graph),
            authors_count=count_authors(graph),
            variants_count=count_variants(graph),
            namespaces_count=count_namespaces(graph),
            # Errors
            undefined_namespaces=get_undefined_namespaces(graph),
            undefined_annotations=get_undefined_annotations(graph),
            namespaces_with_incorrect_names=get_namespaces_with_incorrect_names(graph),
            unused_namespaces=get_unused_namespaces(graph),
            unused_annotations=get_unused_annotations(graph),
            unused_list_annotation_values=get_unused_list_annotation_values(graph),
            naked_names=get_naked_names(graph),
            error_count=count_error_types(graph),
            error_groups=get_most_common_errors(graph),
            syntax_errors=get_syntax_errors(graph),
            # Node pairs
            regulatory_pairs=get_regulatory_pairs(graph),
            chaotic_pairs=get_chaotic_pairs(graph),
            dampened_pairs=get_dampened_pairs(graph),
            contradictory_pairs=get_contradiction_summary(graph),
            separate_unstable_correlation_triples=get_separate_unstable_correlation_triples(graph),
            mutually_unstable_correlation_triples=get_mutually_unstable_correlation_triples(graph),
            jens_unstable=get_jens_unstable(graph),
            increase_mismatch_triplets=get_increase_mismatch_triplets(graph),
            decrease_mismatch_triplets=get_decrease_mismatch_triplets(graph),
            # Bibliometrics
            citation_years=get_citation_years(graph),
            confidence_count=count_confidences(graph),
            # Random
            hub_data=_count_top_hubs(graph),
            disease_data=_count_top_diseases(graph),
        )


def _count_top_hubs(graph, n: int = 15):
    return collections.Counter({
        (
            node.name or node.identifier
            if NAME in node or IDENTIFIER in node else
            str(node)
        ): degree
        for node, degree in get_top_hubs(graph, n=n)
        if isinstance(node, BaseAbundance)
    })


def _count_top_diseases(graph, n=15):
    return collections.Counter({
        (
            node.name or node.identifier
            if NAME in node or IDENTIFIER in node else
            str(node)
        ): count
        for node, count in get_top_pathologies(graph, n=n)
        if isinstance(node, BaseAbundance)
    })
