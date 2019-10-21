# -*- coding: utf-8 -*-

"""Utilities to assemble a BEL graph as bipartite graph of nodes and reified edges."""

import logging
from abc import ABC, abstractmethod
from collections import defaultdict
from itertools import product, starmap
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx

from pybel import BELGraph
from pybel.constants import (
    ACTIVITY, CAUSAL_DECREASE_RELATIONS, CAUSAL_INCREASE_RELATIONS,
    CAUSAL_RELATIONS, DEGRADATION, HAS_COMPONENT, HAS_VARIANT, IS_A, MODIFIER,
    OBJECT, REGULATES, TRANSCRIBED_TO, TRANSLATED_TO, TWO_WAY_RELATIONS,
    BIOMARKER_FOR, CAUSES_NO_CHANGE, HAS_MEMBER
)
from pybel.dsl import BaseEntity, BiologicalProcess, ComplexAbundance, fragment, gene, mirna, pmod, protein, rna
from pybel.struct.filters.edge_predicate_builders import build_relation_predicate
from pybel.struct.filters import filter_edges
from pybel.struct.mutation.collapse.protein_rna_origins import _build_collapse_to_gene_dict
from pybel.typing import EdgeData

__all__ = [
    'reify_bel_graph',
]

logger = logging.getLogger(__name__)

REIF_SUBJECT = 'subject'
REIF_OBJECT = 'object'

ACTIVATES = 'activates'
COMPLEX = 'hasComponent'
CORRELATES = 'correlatesTo'
FRAGMENTS = 'fragments'
INCREASES_ABUNDANCE = 'abundance'
DEGRADATES = 'degradates'
PROMOTES_TRANSLATION = 'translates'
TRANSCRIBES_TO = 'transcribesTo'
TRANSLATES_TO = 'translatesTo'

# Predicates for PTMs

PHOSPHORYLATES = 'phosphorylates'
HYDROXYLATES = 'hydroxylates'
UBIQUITINATES = 'ubiquitinates'
ACETYLATION = "acetylates"
ADP_RIBOSYLATION = "ADP - ribosylates"
FARNESYLATION = "farnesylates"
GERANYLGERANYLATION = "geranylgeranylates"
GLYCOSYLATION = "glycosylates"
ISGYLATION = "ISGylates"
METHYLATION = "methylates"
MYRISTOYLATION = "myristoylates"
NEDDYLATION = "neddylates"
NGLYCO = "glycosylates"
OGLYCO = "glycosylates"
NITROSYLATION = "Nitrosylates"
PALMITOYLATION = "palmitoylates"
SULPHATION = "sulphates"
SUMOYLATION = "SUMOylates"


class ReifiedConverter(ABC):
    """Base class for BEL -> Reified edges graph conversion."""

    target_relation = ...

    @staticmethod
    @abstractmethod
    def predicate(u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        """Test if a BEL edge corresponds to the converter."""

    @staticmethod
    def is_causal_increase(edge_data: EdgeData) -> bool:
        """Check if the relation is ``->``, ``=>``, or ``reg``."""
        return "relation" in edge_data and edge_data['relation'] in CAUSAL_INCREASE_RELATIONS | {REGULATES}

    @staticmethod
    def is_causal_decrease(edge_data: EdgeData) -> bool:
        """Check if the relation is ``-|``, ``=|``, or ``reg``."""
        return "relation" in edge_data and edge_data['relation'] in CAUSAL_DECREASE_RELATIONS | {REGULATES}

    @classmethod
    def convert(
        cls,
        u: BaseEntity,
        v: BaseEntity,
        key: str,
        edge_data: EdgeData
    ) -> Optional[Tuple[BaseEntity, str, bool, bool, BaseEntity]]:
        """Convert a BEL edge to a reified edge.

        Increase and decrease relations have same label, but different sign (positive and negative respectively).
        """
        return (
            u,
            cls.target_relation,
            cls.is_causal_increase(edge_data),
            cls.is_causal_decrease(edge_data),
            v,
        )


class NonCausalConverter(ReifiedConverter):
    """Converts BEL statements of the form A pos B, A neg B and A cor B."""

    target_relation = CORRELATES

    @classmethod
    def predicate(
            cls,
            u: BaseEntity,
            v: BaseEntity,
            key: str,
            edge_data: EdgeData
    ) -> bool:
        return ("relation" in edge_data and
                edge_data['relation'] in TWO_WAY_RELATIONS)


class PTMConverter(ReifiedConverter):
    """Converts BEL statements of the form A X p(B, pmod(*))."""

    # Full list of pmod may be found at
    # http://openbel.org/language/version_2.0/bel_specification_version_2.0.html#_modification_types_provided_in_default_bel_namespace
    synonyms = ...

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            "relation" in edge_data
            and edge_data['relation'] in CAUSAL_RELATIONS
            and "variants" in v
            and any(
                var_['identifier']['name'] in cls.synonyms
                for var_ in v["variants"]
                if isinstance(var_, pmod)
            )
        )


class AbundanceConverter(ReifiedConverter):
    """Convert BEL statements about increasing abundance.

    Statements are of the form ``A B C``, where B in {CAUSAL RELATIONS}
    and A and C don't fall in another special case (pmod, act, ...).
    """

    target_relation = INCREASES_ABUNDANCE

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return "relation" in edge_data and edge_data['relation'] in CAUSAL_RELATIONS


class AcetylationConverter(PTMConverter):
    """Converts BEL statements of the form ``A B p(C, pmod(Ac))``, or synonyms."""

    synonyms = ["Ac", "acetylation"]
    target_relation = ACETYLATION


class ADPRibConverter(PTMConverter):
    """Converts BEL statements of the form ``A B p(C, pmod(ADPRib))``, or synonyms."""

    synonyms = [
        "ADPRib", "ADP - ribosylation",
        "ADP - rybosylation",
        "adenosine diphosphoribosyl",
    ]
    target_relation = ADP_RIBOSYLATION


class ActivationConverter(ReifiedConverter):
    """Converts BEL statements of the form A B act(C)."""

    target_relation = "activates"

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            "relation" in edge_data
            and edge_data['relation'] in CAUSAL_INCREASE_RELATIONS
            and edge_data.get(OBJECT)
            and edge_data.get(OBJECT).get(MODIFIER) == ACTIVITY
        )


class BioProcessConverter(ReifiedConverter):
    """Converts BEL statements of the form A B bp(C)."""

    target_relation = "regulatesBioProcess"

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            "relation" in edge_data and
            edge_data['relation'] in CAUSAL_RELATIONS and
            isinstance(v, BiologicalProcess)
        )


class ComplexConverter(ReifiedConverter):
    """Converts BEL statements of the form complex(A [, **B]])."""

    target_relation = "hasComponent"

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            isinstance(u, ComplexAbundance)
            and "relation" in edge_data
            and edge_data['relation'] in HAS_COMPONENT
        )


class DegradationConverter(ReifiedConverter):
    """Converts BEL statements of the form A B deg(C), when B in {CAUSAL RELATIONS}."""

    target_relation = DEGRADATES

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            "relation" in edge_data
            and edge_data['relation'] in CAUSAL_INCREASE_RELATIONS
            and edge_data.get(OBJECT)
            and edge_data.get(OBJECT).get(MODIFIER) == DEGRADATION
        )


class FarnesylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Farn)) or synonyms."""

    synonyms = ["Farn", "farnesylation"]
    target_relation = FARNESYLATION


class GeranylgeranylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Gerger)) or synonyms."""

    synonyms = ["Gerger", "geranylgeranylation"]
    target_relation = GERANYLGERANYLATION


class HasVariantConverter(ReifiedConverter):
    """Identifies edges of the form A hasVariant B.

    Do not convert them to reified edges.
    """

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return "relation" in edge_data and edge_data['relation'] == HAS_VARIANT

    @classmethod
    def convert(
        cls,
        u: BaseEntity,
        v: BaseEntity,
        key: str,
        edge_data: EdgeData,
    ) -> Optional[Tuple[BaseEntity, str, bool, bool, BaseEntity]]:
        return None


class IsAConverter(ReifiedConverter):
    """Identifies edges of the form A isA B. Do not convert them to
    reified edges."""

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return ("relation" in edge_data and
                edge_data['relation'] == IS_A)

    @classmethod
    def convert(cls,
                u: BaseEntity,
                v: BaseEntity,
                key: str,
                edge_data: EdgeData
                ) -> Optional[Tuple[BaseEntity, str, bool, bool, BaseEntity]]:
        return None


class NonReifiableRelationshipsConverter(ReifiedConverter):
    """Identifies edges that won`t be reified."""

    do_not_reify = [BIOMARKER_FOR, CAUSES_NO_CHANGE, HAS_MEMBER]

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return ("relation" in edge_data and
                edge_data['relation'] in cls.do_not_reify)

    @classmethod
    def convert(cls,
                u: BaseEntity,
                v: BaseEntity,
                key: str,
                edge_data: EdgeData
                ) -> Optional[Tuple[BaseEntity, str, bool, bool, BaseEntity]]:
        return None


class FragmentationConverter(ReifiedConverter):
    """Converts BEL statements of the form A B p(C, frag(*))."""

    target_relation = FRAGMENTS

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            "relation" in edge_data
            and edge_data['relation'] in CAUSAL_RELATIONS
            and "variants" in v
            and any(
                isinstance(var_, fragment)
                for var_ in v["variants"]
            )
        )


class GlycosylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Glyco)) or synonyms."""

    synonyms = [
        "Glyco", "glycosylation", "NGlyco", "N - linked glycosylation",
        "OGlyco", "O - linked glycosylation"
    ]
    target_relation = GLYCOSYLATION


class HydroxylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Hy)) or synonyms."""

    synonyms = ["Hy", "hydroxylation"]
    target_relation = HYDROXYLATES


class ISGylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(ISG)) or synonyms."""

    synonyms = ["ISG", "ISGylation", "ISG15 - protein conjugation"]
    target_relation = ISGYLATION


class MethylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Me)) or synonyms."""

    synonyms = [
        "Me", "methylation", "Me1", "monomethylation", "mono - methylation",
        "Me2", "dimethylation", "di - methylation", "Me3", "trimethylation",
        "tri - methylation"
    ]
    target_relation = METHYLATION


class MyristoylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Myr)) or synonyms."""

    synonyms = ["Myr", "myristoylation"]
    target_relation = MYRISTOYLATION


class NeddylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Nedd)) or synonyms."""

    synonyms = ["Nedd", "neddylation"]
    target_relation = NEDDYLATION


class NitrosylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(NO)) or synonyms."""

    synonyms = ["NO", "Nitrosylation"]
    target_relation = NITROSYLATION


class PalmitoylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Palm)) or synonyms."""

    synonyms = ["Palm", "palmitoylation"]
    target_relation = PALMITOYLATION


class PhosphorylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Ph)) or synonyms."""

    synonyms = ["Ph", "phosphorylation"]
    target_relation = PHOSPHORYLATES


class PromotesTranslationConverter(ReifiedConverter):
    """Converts BEL statements of the form A X r(B)."""

    target_relation = PROMOTES_TRANSLATION

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            "relation" in edge_data
            and edge_data['relation'] in CAUSAL_RELATIONS
            and isinstance(v, rna)
        )


class SulphationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Sulf)) or synonyms."""

    synonyms = [
        "Sulf", "sulfation", "sulphation", "sulfur addition",
        "sulphur addition", "sulfonation", "sulphonation"
    ]
    target_relation = SULPHATION


class SUMOylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(SUMO)) or synonyms."""

    synonyms = ["Sumo", "SUMOylation"]
    target_relation = SUMOYLATION


class TranscriptionConverter(ReifiedConverter):
    """Converts BEL statements of the form g(A) :> r(C)."""

    target_relation = PROMOTES_TRANSLATION

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            "relation" in edge_data
            and edge_data['relation'] in TRANSCRIBED_TO
            and isinstance(u, gene)
            and isinstance(v, rna)
        )


class TranslationConverter(ReifiedConverter):
    """Converts BEL statements of the form r(A) >> p(C)."""

    target_relation = PROMOTES_TRANSLATION

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            "relation" in edge_data
            and edge_data['relation'] in TRANSLATED_TO
            and isinstance(u, rna)
            and isinstance(v, protein)
        )


class UbiquitinationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Ub)) or synonyms."""

    synonyms = [
        "Ub", "ubiquitination", "ubiquitinylation", "ubiquitylation", "UbK48",
        "Lysine 48 - linked polyubiquitination", "UbK63", "UbMono", "UbPoly",
        "Lysine 63 - linked polyubiquitination", "monoubiquitination"
    ]
    target_relation = UBIQUITINATES


def reify_edge(
    u: BaseEntity,
    v: BaseEntity,
    key: str,
    edge_data: EdgeData,
) -> Optional[Tuple[BaseEntity, str, bool, bool, BaseEntity]]:
    """Reify a given edge, if possible."""
    converters = [
        TranslationConverter,
        TranscriptionConverter,
        ComplexConverter,
        FragmentationConverter,
        ADPRibConverter,
        FarnesylationConverter,
        GeranylgeranylationConverter,
        GlycosylationConverter,
        ISGylationConverter,
        MethylationConverter,
        MyristoylationConverter,
        NeddylationConverter,
        NitrosylationConverter,
        PalmitoylationConverter,
        PhosphorylationConverter,
        UbiquitinationConverter,
        SulphationConverter,
        SUMOylationConverter,
        HydroxylationConverter,
        ActivationConverter,
        DegradationConverter,
        BioProcessConverter,
        PromotesTranslationConverter,
        AbundanceConverter,
        HasVariantConverter,
        NonCausalConverter,
        IsAConverter,
        NonReifiableRelationshipsConverter
    ]
    for converter in converters:
        if converter.predicate(u, v, key, edge_data):
            return converter.convert(u, v, key, edge_data)

    logger.warning(f"No converter found for {u}, {v}")
    logger.warning(f"  with edge data {edge_data}")


def _get_entity_collapser(bel_graph: BELGraph, collapse: str) -> Dict[BaseEntity, Set[BaseEntity]]:
    """Creates a mapping to collapsed entities from a BEL graph.

    :param bel_graph: The BEL graph.
    :param collapse: How to collapse the entities: ``variants`` to remove the variants, ``genes`` to collapse protein,
     RNA, and miRNA nodes to their corresponding gene nodes or None to maintain the graph as is.
    :return: A mapping from the original entity to the collapsed entity.
    """
    if collapse == 'variants':
        has_variant_predicate = build_relation_predicate(HAS_VARIANT)
        edges = filter_edges(bel_graph, has_variant_predicate)
        collapse_dict = defaultdict(set)
        for orig, vari, _ in edges:
            collapse_dict[vari].add(orig)
        return collapse_dict
    elif collapse == 'genes':
        collapse_dict = _build_collapse_to_gene_dict(bel_graph)
        return collapse_dict
    else:
        return dict()


def _collapse_edges(
        edges: List[Tuple[BaseEntity, str, bool, bool, BaseEntity]],
        collapser: Dict[BaseEntity, Set[BaseEntity]]
) -> Set[Tuple[BaseEntity, str, bool, bool, BaseEntity]]:
    """Collapses the entities in a list with edge components and remove duplicate relationships.

    :param edges: A list of edges, which are represented in a tuple containing: the origin node, the predicate,
    positive causal relationship flag, negative causal relationship flag and the target node.
    :param collapser: A mapping directing the collapse.
    :return: A set of edges with collapsed entities.
    """
    collapsed_edges = set()
    for source_node, reif_edge_label, positive, negative, target_node in edges:
        if {source_node, target_node}.intersection(collapser):
            source_set = collapser[source_node] if source_node in collapser else [source_node]
            target_set = collapser[target_node] if target_node in collapser else [target_node]
            collapsed_edges = collapsed_edges.union(
                product(source_set, [reif_edge_label], [positive], [negative], target_set)
            )
        else:
            collapsed_edges.add((source_node, reif_edge_label, positive, negative, target_node))
    return collapsed_edges


def reify_bel_graph(bel_graph: BELGraph, collapse: str = None) -> nx.DiGraph:
    """Generate a new graph with reified edges."""
    reified_graph = nx.DiGraph()

    reified_edges = starmap(reify_edge, bel_graph.edges(keys=True, data=True))
    reified_edges = list(filter(None, reified_edges))

    collapser = _get_entity_collapser(bel_graph, collapse)
    collapsed_edges = _collapse_edges(reified_edges, collapser)

    for i, (source, reif_edge_label, positive, negative, target) in enumerate(collapsed_edges):
        reified_graph.add_node(i, label=reif_edge_label, causal=(positive, negative))
        reified_graph.add_edge(source, i, label=REIF_SUBJECT)
        reified_graph.add_edge(i, target, label=REIF_OBJECT)

    return reified_graph
