# -*- coding: utf-8 -*-

"""Utilities to assemble a BEL graph as bipartite graph of nodes and reified edges."""

import logging
import warnings
from abc import ABC, abstractmethod
from itertools import starmap
from typing import Optional, Tuple

import networkx as nx

import pybel.constants as pbc
from pybel import BELGraph
from pybel.dsl import BaseEntity, CentralDogma, ComplexAbundance, ProteinModification, fragment, gene, protein, rna
from pybel.typing import EdgeData

__all__ = [
    'reify_bel_graph',
]

REIF_SUBJECT = 'subject'
REIF_OBJECT = 'object'

ACTIVATES = 'activates'
COMPLEX = 'hasComponent'
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
        return pbc.RELATION in edge_data and edge_data[pbc.RELATION] in pbc.CAUSAL_INCREASE_RELATIONS | {pbc.REGULATES}

    @staticmethod
    def is_causal_decrease(edge_data: EdgeData) -> bool:
        """Check if the relation is ``-|``, ``=|``, or ``reg``."""
        return pbc.RELATION in edge_data and edge_data[pbc.RELATION] in pbc.CAUSAL_DECREASE_RELATIONS | {pbc.REGULATES}

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


class PTMConverter(ReifiedConverter):
    """Converts BEL statements of the form A X p(B, pmod(*))."""

    synonyms = ...

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            edge_data[pbc.RELATION] in pbc.CAUSAL_RELATIONS
            and isinstance(v, CentralDogma)
            and v.variants
            and any(
                variant.entity.name in cls.synonyms
                for variant in v.variants
                if isinstance(variant, ProteinModification)
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
        return pbc.RELATION in edge_data and edge_data[pbc.RELATION] in pbc.CAUSAL_RELATIONS


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
            pbc.RELATION in edge_data
            and edge_data[pbc.RELATION] in pbc.CAUSAL_INCREASE_RELATIONS
            and edge_data.get(pbc.TARGET_MODIFIER)
            and edge_data.get(pbc.TARGET_MODIFIER).get(pbc.MODIFIER) == pbc.ACTIVITY
        )


class ComplexConverter(ReifiedConverter):
    """Converts BEL statements of the form complex(A [, **B]])."""

    target_relation = "partOf"

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            isinstance(v, ComplexAbundance)
            and edge_data[pbc.RELATION] == pbc.PART_OF
        )


class DegradationConverter(ReifiedConverter):
    """Converts BEL statements of the form A B act(C), when B in {CAUSAL RELATIONS}."""

    target_relation = DEGRADATES

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            pbc.RELATION in edge_data
            and edge_data[pbc.RELATION] in pbc.CAUSAL_INCREASE_RELATIONS
            and edge_data.get(pbc.TARGET_MODIFIER)
            and edge_data.get(pbc.TARGET_MODIFIER).get(pbc.MODIFIER) == pbc.DEGRADATION
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
        return pbc.RELATION in edge_data and edge_data[pbc.RELATION] == pbc.HAS_VARIANT

    @classmethod
    def convert(
        cls,
        u: BaseEntity,
        v: BaseEntity,
        key: str,
        edge_data: EdgeData,
    ) -> Optional[Tuple[BaseEntity, str, bool, bool, BaseEntity]]:
        return None


class FragmentationConverter(ReifiedConverter):
    """Converts BEL statements of the form A B p(C, pmod(Hy))."""

    target_relation = FRAGMENTS

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            pbc.RELATION in edge_data
            and edge_data[pbc.RELATION] in pbc.CAUSAL_RELATIONS
            and pbc.VARIANTS in v
            and any(
                isinstance(var_, fragment)
                for var_ in v[pbc.VARIANTS]
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
            pbc.RELATION in edge_data
            and edge_data[pbc.RELATION] in pbc.CAUSAL_RELATIONS
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
            pbc.RELATION in edge_data
            and edge_data[pbc.RELATION] in pbc.TRANSCRIBED_TO
            and isinstance(u, gene)
            and isinstance(v, rna)
        )


class TranslationConverter(ReifiedConverter):
    """Converts BEL statements of the form r(A) >> p(C)."""

    target_relation = PROMOTES_TRANSLATION

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity, key: str, edge_data: EdgeData) -> bool:
        return (
            pbc.RELATION in edge_data
            and edge_data[pbc.RELATION] in pbc.TRANSLATED_TO
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
        PromotesTranslationConverter,
        AbundanceConverter,
        HasVariantConverter
    ]
    for converter in converters:
        if converter.predicate(u, v, key, edge_data):
            return converter.convert(u, v, key, edge_data)

    logging.warning(f"No converter found for {u}, {v}")
    logging.warning(f"  with edge data {edge_data}")


def reify_bel_graph(bel_graph: BELGraph) -> nx.DiGraph:
    """Generate a new graph with reified edges."""
    warnings.warn(
        'reify_bel_graph does not pass unit tests. Probably does not work.'
        'Use with caution or email author Mauricio de ')

    reified_graph = nx.DiGraph()

    reified_edges = starmap(reify_edge, bel_graph.edges(keys=True, data=True))
    reified_edges = filter(None, reified_edges)

    for i, (source, reif_edge_label, positive, negative, target) in enumerate(reified_edges):
        reified_graph.add_node(i, label=reif_edge_label, causal=(positive, negative))
        reified_graph.add_edge(source, i, label=REIF_SUBJECT)
        reified_graph.add_edge(target, i, label=REIF_OBJECT)

    return reified_graph
