# -*- coding: utf-8 -*-

"""Utilities to assemble a BEL graph as bipartite graph of nodes and
reified edges."""


import logging
import unittest
from abc import ABC, abstractmethod
from itertools import count
from typing import Dict, Optional, Tuple

import networkx as nx

from pybel import BELGraph
from pybel.constants import (
    ACTIVITY, CAUSAL_DECREASE_RELATIONS, CAUSAL_INCREASE_RELATIONS,
    CAUSAL_RELATIONS, DEGRADATION, HAS_COMPONENT, HAS_VARIANT, MODIFIER,
    OBJECT, REGULATES, TRANSCRIBED_TO, TRANSLATED_TO
)
from pybel.dsl import (
    abundance, activity, BaseEntity, ComplexAbundance, degradation, fragment,
    gene, pmod, protein, rna
)
from pybel.testing.utils import n

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
    def predicate(u: BaseEntity, v: BaseEntity,
                  key: str, edge_data: Dict) -> bool:
        """Test if a BEL edge corresponds to the converter."""

    @staticmethod
    def is_causal_increase(edge_data: Dict) -> bool:
        """Checks if the relation is ->, => or reg"""

        return ("relation" in edge_data and
                edge_data['relation'] in CAUSAL_INCREASE_RELATIONS | {REGULATES})

    @staticmethod
    def is_causal_decrease(edge_data: Dict) -> bool:
        """Checks if the relation is -|, =| or reg"""

        return ("relation" in edge_data and
                edge_data['relation'] in CAUSAL_DECREASE_RELATIONS | {REGULATES})

    @classmethod
    def convert(cls,
                u: BaseEntity,
                v: BaseEntity,
                key: str,
                edge_data: Dict
                ) -> Optional[Tuple[BaseEntity, str, bool, bool, BaseEntity]]:
        """Convert a BEL edge to a reified edge. Increase and decrease
        relations have same label, but different sign (positive and negative
        respectively)."""

        return (u,
                cls.target_relation,
                cls.is_causal_increase(edge_data),
                cls.is_causal_decrease(edge_data),
                v)


class PTMConverter(ReifiedConverter):
    """Converts BEL statements of the form A X p(B, pmod(*))."""

    synonyms = ...

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity,
                  key: str, edge_data: Dict) -> bool:
        return ("relation" in edge_data and
                edge_data['relation'] in CAUSAL_RELATIONS and
                "variants" in v and
                any([var_['identifier']['name'] in cls.synonyms
                     for var_ in v["variants"]
                     if isinstance(var_, pmod)]))


class AbundanceConverter(ReifiedConverter):
    """Converts BEL statements of the form A B C, where B in {CAUSAL RELATIONS}
    and A and C don't fall in another special case (pmod, act, ...).
    """

    target_relation = INCREASES_ABUNDANCE

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity,
                  key: str, edge_data: Dict) -> bool:
        return ("relation" in edge_data and
                edge_data['relation'] in CAUSAL_RELATIONS)


class AcetylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Ac)) or
    synonyms."""

    synonyms = ["Ac", "acetylation"]
    target_relation = ACETYLATION


class ADPRibConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(ADPRib)) or
    synonyms."""

    synonyms = [
        "ADPRib", "ADP - ribosylation", "ADP - rybosylation",
        "adenosine diphosphoribosyl"
    ]
    target_relation = ADP_RIBOSYLATION


class ActivationConverter(ReifiedConverter):
    """Converts BEL statements of the form A B act(C)."""

    target_relation = "activates"

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity,
                  key: str, edge_data: Dict) -> bool:
        return (
                "relation" in edge_data and
                edge_data['relation'] in CAUSAL_INCREASE_RELATIONS and
                edge_data.get(OBJECT) and
                edge_data.get(OBJECT).get(MODIFIER) == ACTIVITY
        )


class ComplexConverter(ReifiedConverter):
    """Converts BEL statements of the form complex(A [, **B]])."""

    target_relation = "hasComponent"

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity,
                  key: str, edge_data: Dict) -> bool:
        return (
            isinstance(u, ComplexAbundance) and
            "relation" in edge_data and
            edge_data['relation'] in HAS_COMPONENT
        )


class DegradationConverter(ReifiedConverter):
    """Converts BEL statements of the form A B act(C), when B in
    {CAUSAL RELATIONS}."""

    target_relation = DEGRADATES

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity,
                  key: str, edge_data: Dict) -> bool:
        return (
                "relation" in edge_data and
                edge_data['relation'] in CAUSAL_INCREASE_RELATIONS and
                edge_data.get(OBJECT) and
                edge_data.get(OBJECT).get(MODIFIER) == DEGRADATION
        )


class FarnesylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Farn)) or
    synonyms."""

    synonyms = ["Farn", "farnesylation"]
    target_relation = FARNESYLATION


class GeranylgeranylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Gerger)) or
    synonyms."""

    synonyms = ["Gerger", "geranylgeranylation"]
    target_relation = GERANYLGERANYLATION


class HasVariantConverter(ReifiedConverter):
    """Identifies edges of the form A hasvariant B. Do not convert them to
    reified edges."""

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity,
                  key: str, edge_data: Dict) -> bool:
        return ("relation" in edge_data and
                edge_data['relation'] == HAS_VARIANT)

    @classmethod
    def convert(cls,
                u: BaseEntity,
                v: BaseEntity,
                key: str,
                edge_data: Dict
                ) -> Optional[Tuple[BaseEntity, str, bool, bool, BaseEntity]]:
        return None


class FragmentationConverter(ReifiedConverter):
    """Converts BEL statements of the form A B p(C, pmod(Hy))."""

    target_relation = FRAGMENTS

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity,
                  key: str, edge_data: Dict) -> bool:
        return (
                "relation" in edge_data and
                edge_data['relation'] in CAUSAL_RELATIONS and
                "variants" in v and
                any([isinstance(var_, fragment) for var_ in v["variants"]])
        )


class GlycosylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Glyco)) or
    synonyms."""

    synonyms = [
        "Glyco", "glycosylation", "NGlyco", "N - linked glycosylation",
        "OGlyco", "O - linked glycosylation"
    ]
    target_relation = GLYCOSYLATION


class HydroxylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Hy)) or
    synonyms."""

    synonyms = ["Hy", "hydroxylation"]
    target_relation = HYDROXYLATES


class ISGylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(ISG)) or
    synonyms."""

    synonyms = ["ISG", "ISGylation", "ISG15 - protein conjugation"]
    target_relation = ISGYLATION


class MethylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Me)) or
    synonyms."""

    synonyms = [
        "Me", "methylation", "Me1", "monomethylation", "mono - methylation",
        "Me2", "dimethylation", "di - methylation", "Me3", "trimethylation",
        "tri - methylation"
    ]
    target_relation = METHYLATION


class MyristoylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Myr)) or
    synonyms."""

    synonyms = ["Myr", "myristoylation"]
    target_relation = MYRISTOYLATION


class NeddylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Nedd)) or
    synonyms."""

    synonyms = ["Nedd", "neddylation"]
    target_relation = NEDDYLATION


class NitrosylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(NO)) or
    synonyms."""

    synonyms = ["NO", "Nitrosylation"]
    target_relation = NITROSYLATION


class PalmitoylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Palm)) or
    synonyms."""

    synonyms = ["Palm", "palmitoylation"]
    target_relation = PALMITOYLATION


class PhosphorylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Ph)) or
    synonyms."""

    synonyms = ["Ph", "phosphorylation"]
    target_relation = PHOSPHORYLATES


class PromotesTranslationConverter(ReifiedConverter):
    """Converts BEL statements of the form A X r(B)."""

    target_relation = PROMOTES_TRANSLATION

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity,
                  key: str, edge_data: Dict) -> bool:
        return ("relation" in edge_data and
                edge_data['relation'] in CAUSAL_RELATIONS and
                isinstance(v, rna))


class SulphationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Sulf)) or
    synonyms."""

    synonyms = [
        "Sulf", "sulfation", "sulphation", "sulfur addition",
        "sulphur addition", "sulfonation", "sulphonation"
    ]
    target_relation = SULPHATION


class SUMOylationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(SUMO)) or
    synonyms."""

    synonyms = ["Sumo", "SUMOylation"]
    target_relation = SUMOYLATION


class TranscriptionConverter(ReifiedConverter):
    """Converts BEL statements of the form g(A) :> r(C)."""

    target_relation = PROMOTES_TRANSLATION

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity,
                  key: str, edge_data: Dict) -> bool:
        return ("relation" in edge_data and
                edge_data['relation'] in TRANSCRIBED_TO and
                isinstance(u, gene) and
                isinstance(v, rna))


class TranslationConverter(ReifiedConverter):
    """Converts BEL statements of the form r(A) >> p(C)."""

    target_relation = PROMOTES_TRANSLATION

    @classmethod
    def predicate(cls, u: BaseEntity, v: BaseEntity,
                  key: str, edge_data: Dict) -> bool:
        return ("relation" in edge_data and
                edge_data['relation'] in TRANSLATED_TO and
                isinstance(u, rna) and
                isinstance(v, protein)
                )


class UbiquitinationConverter(PTMConverter):
    """Converts BEL statements of the form A B p(C, pmod(Ub)) or
    synonyms."""

    synonyms = [
        "Ub", "ubiquitination", "ubiquitinylation", "ubiquitylation", "UbK48",
        "Lysine 48 - linked polyubiquitination", "UbK63", "UbMono", "UbPoly",
        "Lysine 63 - linked polyubiquitination", "monoubiquitination"
    ]
    target_relation = UBIQUITINATES


def reify_edge(u: BaseEntity,
               v: BaseEntity,
               key: str,
               edge_data: Dict
               ) -> Optional[Tuple[BaseEntity, str, bool, bool, BaseEntity]]:
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

    # No converter found
    return None


def reify_bel_graph(bel_graph: BELGraph) -> nx.DiGraph:
    """Generate a new graph with reified edges."""
    reified_graph = nx.DiGraph()
    gen = count()

    for edge in bel_graph.edges(keys=True):
        (u, v, key) = edge
        data = bel_graph[u][v][key]

        reified_edge = reify_edge(u, v, key, data)
        if reified_edge:
            new_u, reif_edge_label, positive, negative, new_v = reified_edge
            reif_edge_num = next(gen)
            causal = (positive, negative)
            reified_graph.add_node(
                reif_edge_num, label=reif_edge_label, causal=causal
            )
            reified_graph.add_edge(new_u, reif_edge_num, label=REIF_SUBJECT)
            reified_graph.add_edge(new_v, reif_edge_num, label=REIF_OBJECT)

    return reified_graph


cdk5 = protein('HGNC', 'CDK5', 'HGNC:1774')
gsk3b = protein('HGNC', 'GSK3B', 'HGNC:4617')
p_tau = protein('HGNC', 'MAPT', 'HGNC:6893', variants=pmod('Ph'))

# act(p(HGNC:FAS), ma(cat)) increases act(p(HGNC:CASP8), ma(cat))
fas = protein('HGNC', 'FAS', 'HGNC:11920')
casp8 = protein('HGNC', 'CASP8', 'HGNC:1509')

# a(CHEBI:oxaliplatin) increases a(MESHC:"Reactive Oxygen Species")
oxaliplatin = abundance('CHEBI', 'oxaliplatin', 'CHEBI:31941')
reactive_o_species = abundance('MESHC', 'Reactive Oxygen Species', 'D017382')

# p(HGNC:MYC) decreases r(HGNC:CCNB1)


class TestAssembleReifiedGraph(unittest.TestCase):
    """Test assembly of reified graphs."""

    help_causal_increases = (True, False)
    help_causal_decreases = (False, True)
    help_not_causal = (False, False)
    help_causal_regulates = (True, True)

    def help_test_graphs_equal(self,
                               expected: nx.DiGraph,
                               actual: nx.DiGraph
                               ) -> None:
        """Test that two DiGraphs are equal."""
        self.assertIsNotNone(actual)
        self.assertIsInstance(actual, nx.DiGraph)
        self.assertEqual(expected.number_of_nodes(), actual.number_of_nodes())
        self.assertEqual(expected.number_of_edges(), actual.number_of_edges())

        for node in expected:
            self.assertIn(node, actual)

        actual_edges_list = [(u_, actual.nodes[v_]['label'], actual.nodes[v_]['causal'])
                             for u_, v_ in actual.edges]

        for u, v in expected.edges():
            self.assertIn((u, expected.nodes[v]['label'], expected.nodes[v]['causal']), actual_edges_list)

    def test_convert_dephosphorylates(self):
        """Test the conversion of a BEL statement like
        ``act(p(X)) -| p(Y, pmod(Ph))."""
        bel_graph = BELGraph()
        bel_graph.add_directly_decreases(
            cdk5,
            p_tau,
            evidence=n(),
            citation=n(),
            subject_modifier=activity('kin'),
        )

        r_edge = 0
        expected_reified_graph = \
            TestAssembleReifiedGraph.help_make_simple_expected_graph(
                cdk5,
                p_tau,
                PHOSPHORYLATES,
                r_edge,
                self.help_causal_decreases
            )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_phosphorylates(self):
        """Test the conversion of a BEL statement like
        ``act(p(X)) -> p(Y, pmod(Ph))."""
        bel_graph = BELGraph()
        bel_graph.add_directly_increases(
            cdk5,
            p_tau,
            evidence=n(),
            citation=n(),
            subject_modifier=activity('kin'),
        )

        r_edge = 0
        expected_reified_graph = \
            TestAssembleReifiedGraph.help_make_simple_expected_graph(
                cdk5,
                p_tau,
                PHOSPHORYLATES,
                r_edge,
                self.help_causal_increases
            )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_two_phosphorylates(self):
        """Test that two phosphorylations of the same object get
        different reified nodes."""
        bel_graph = BELGraph()
        for kinase in (cdk5, gsk3b):
            bel_graph.add_directly_increases(
                kinase,
                p_tau,
                evidence=n(),
                citation=n(),
                subject_modifier=activity('kin'),
            )

        re1, re2 = 0, 1
        expected_reified_graph = \
            TestAssembleReifiedGraph.help_make_simple_expected_graph(
                cdk5,
                p_tau,
                PHOSPHORYLATES,
                re1,
                self.help_causal_increases
            )
        expected_reified_graph.add_node(
            re2, label='phosphorylates', causal=self.help_causal_increases
        )
        expected_reified_graph.add_edge(gsk3b, re2, label=REIF_SUBJECT)
        expected_reified_graph.add_edge(p_tau, re2, label=REIF_OBJECT)

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_activates(self):
        """Test the conversion of a bel statement like p(x) -> act(p(y))."""

        bel_graph = BELGraph()
        bel_graph.add_directly_increases(
            cdk5,
            casp8,
            evidence=n(),
            citation=n(),
            object_modifier=activity('ma')
        )

        expected_reified_graph = \
            TestAssembleReifiedGraph.help_make_simple_expected_graph(
                cdk5,
                casp8,
                ACTIVATES,
                0,
                self.help_causal_increases
            )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_increases_abundance(self):
        """Test the conversion of a bel statement like A -> B, when
        A and B don't fall in any special case (activity, pmod, ...).
        """

        bel_graph = BELGraph()
        bel_graph.add_increases(
            oxaliplatin,
            reactive_o_species,
            evidence='10.1093/jnci/djv394',
            citation='PubMed:26719345'
        )

        expected_reified_graph = \
            TestAssembleReifiedGraph.help_make_simple_expected_graph(
                oxaliplatin,
                reactive_o_species,
                INCREASES_ABUNDANCE,
                0,
                self.help_causal_increases
            )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_degradates(self):
        """Test the conversion of a bel statement like A -> deg(B).
        """

        microglia = abundance('MeSH', 'Microglia', 'MeSH:D017628')
        abeta = abundance('CHEBI', 'amyloid-β', 'CHEBI:64645')

        # a(MESH:Microglia) reg deg(a(CHEBI:"amyloid-beta"))
        bel_graph = BELGraph()
        bel_graph.add_increases(
            microglia,
            abeta,
            evidence='10.1038/s41586-018-0368-8',
            citation='PubMed:30046111',
            object_modifier=degradation()
        )

        expected_reified_graph = \
            TestAssembleReifiedGraph.help_make_simple_expected_graph(
                microglia,
                abeta,
                DEGRADATES,
                0,
                self.help_causal_increases
            )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_fragments(self):
        """Test the conversion of a bel statement like A -> p(B, frag(?)).
        """

        casp3 = abundance('HGNC', 'CASP3', 'MeSH:D017628')
        tau = protein('HGNC', 'MAPT', 'HGNC:6893', variants=fragment())

        # act(p(HGNC:CASP3), ma(pep)) increases p(HGNC:MAPT, frag("?"))
        bel_graph = BELGraph()
        bel_graph.add_increases(
            casp3,
            tau,
            evidence='10.1038/s41586-018-0368-8',
            citation='PubMed:30046111'
        )

        expected_reified_graph = \
            TestAssembleReifiedGraph.help_make_simple_expected_graph(
                casp3,
                tau,
                FRAGMENTS,
                0,
                self.help_causal_increases
            )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_promote_translation(self):
        """Test the conversion of a bel statement like A -> r(B)"""

        # example from Colorectal Cancer Model v2.0.6 @ scai
        # act(p(HGNC:CTNNB1), ma(tscript)) increases r(HGNC:BIRC5)
        ctnnb1 = protein('HGNC', 'CTNNB1', '')
        birc5 = rna('HGNC', 'BIRC5', '')

        # a(MESH:Microglia) reg deg(a(CHEBI:"amyloid-beta"))
        bel_graph = BELGraph()
        bel_graph.add_increases(
            ctnnb1,
            birc5,
            evidence='10.1038/s41586-018-0368-8',
            citation='PMID:18075512',
            subject_modifier=activity('tscript')
        )

        expected_reified_graph = \
            TestAssembleReifiedGraph.help_make_simple_expected_graph(
                ctnnb1,
                birc5,
                PROMOTES_TRANSLATION,
                0,
                self.help_causal_increases
            )
        reified_graph = reify_bel_graph(bel_graph)

        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_increases_abundance_then_phosphorylates(self):
        """Test the conversion of a bel graph containing one increases
        abundance and one phosphorylates relationship."""

        bel_graph = BELGraph()
        bel_graph.add_increases(
            oxaliplatin,
            reactive_o_species,
            evidence='10.1093/jnci/djv394',
            citation='PubMed:26719345'
        )
        bel_graph.add_directly_increases(
            reactive_o_species,
            p_tau,
            evidence=n(),
            citation=n()
        )

        re1, re2 = 1, 0
        expected_reified_graph = \
            TestAssembleReifiedGraph.help_make_simple_expected_graph(
                oxaliplatin,
                reactive_o_species,
                INCREASES_ABUNDANCE,
                re1,
                self.help_causal_increases
            )

        expected_reified_graph.add_node(
            re2, label=PHOSPHORYLATES, causal=self.help_causal_increases
        )
        expected_reified_graph.add_edge(
            reactive_o_species, re2, label=REIF_SUBJECT
        )
        expected_reified_graph.add_edge(
            p_tau, re2, label=REIF_OBJECT
        )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    @staticmethod
    def help_make_simple_expected_graph(u, v, label, edge_num, causal_tup):
        """Creates a simple reified graph (with 3 nodes - subject, predicate
        and object)."""

        expected_reified_graph = nx.DiGraph()
        expected_reified_graph.add_node(edge_num, label=label, causal=causal_tup)
        expected_reified_graph.add_edge(u, edge_num, label=REIF_SUBJECT)
        expected_reified_graph.add_edge(v, edge_num, label=REIF_OBJECT)
        return expected_reified_graph


if __name__ == '__main__':
    unittest.main()
