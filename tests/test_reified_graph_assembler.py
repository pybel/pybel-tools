# -*- coding: utf-8 -*-

"""Tests for the reified assembler."""

import unittest

import networkx as nx
from pybel_tools.assembler.reified_graph.assembler import (
    ACTIVATES, DEGRADATES, FRAGMENTS, INCREASES_ABUNDANCE, PHOSPHORYLATES, PROMOTES_TRANSLATION, REIF_OBJECT,
    REIF_SUBJECT, reify_bel_graph,
)

from pybel import BELGraph
from pybel.dsl import abundance, activity, degradation, fragment, pmod, protein, rna
from pybel.testing.utils import n

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

    def help_test_graphs_equal(
            self,
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

        actual_edges_list = [
            (u_, actual.nodes[v_]['label'], actual.nodes[v_]['causal'])
            for u_, v_ in actual.edges
        ]

        for u, v in expected.edges():
            self.assertIn((u, expected.nodes[v]['label'], expected.nodes[v]['causal']), actual_edges_list)

    def test_convert_dephosphorylates(self):
        """Test the conversion of a BEL statement like ``act(p(X)) -| p(Y, pmod(Ph))."""
        bel_graph = BELGraph()
        bel_graph.add_directly_decreases(
            cdk5,
            p_tau,
            evidence=n(),
            citation=n(),
            subject_modifier=activity('kin'),
        )

        r_edge = 0
        expected_reified_graph = self.help_make_simple_expected_graph(
            cdk5,
            p_tau,
            PHOSPHORYLATES,
            r_edge,
            self.help_causal_decreases,
        )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_phosphorylates(self):
        """Test the conversion of a BEL statement like ``act(p(X)) -> p(Y, pmod(Ph))."""
        bel_graph = BELGraph()
        bel_graph.add_directly_increases(
            cdk5,
            p_tau,
            evidence=n(),
            citation=n(),
            subject_modifier=activity('kin'),
        )

        r_edge = 0
        expected_reified_graph = self.help_make_simple_expected_graph(
            cdk5,
            p_tau,
            PHOSPHORYLATES,
            r_edge,
            self.help_causal_increases
        )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_two_phosphorylates(self):
        """Test that two phosphorylations of the same object get different reified nodes."""
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
        expected_reified_graph = self.help_make_simple_expected_graph(
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

        expected_reified_graph = self.help_make_simple_expected_graph(
            cdk5,
            casp8,
            ACTIVATES,
            0,
            self.help_causal_increases
        )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_increases_abundance(self):
        """Test the conversion of a bel statement like A -> B, when A and B don't fall in any special case (activity, pmod, ...)."""
        bel_graph = BELGraph()
        bel_graph.add_increases(
            oxaliplatin,
            reactive_o_species,
            evidence='10.1093/jnci/djv394',
            citation='PubMed:26719345',
        )

        expected_reified_graph = self.help_make_simple_expected_graph(
            oxaliplatin,
            reactive_o_species,
            INCREASES_ABUNDANCE,
            0,
            self.help_causal_increases,
        )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_degradates(self):
        """Test the conversion of a bel statement like A -> deg(B)."""

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

        expected_reified_graph = self.help_make_simple_expected_graph(
            microglia,
            abeta,
            DEGRADATES,
            0,
            self.help_causal_increases
        )

        reified_graph = reify_bel_graph(bel_graph)
        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_fragments(self):
        """Test the conversion of a bel statement like A -> p(B, frag(?))."""

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

        expected_reified_graph = self.help_make_simple_expected_graph(
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

        expected_reified_graph = self.help_make_simple_expected_graph(
            ctnnb1,
            birc5,
            PROMOTES_TRANSLATION,
            0,
            self.help_causal_increases
        )
        reified_graph = reify_bel_graph(bel_graph)

        self.help_test_graphs_equal(expected_reified_graph, reified_graph)

    def test_convert_increases_abundance_then_phosphorylates(self):
        """Test the conversion of a bel graph containing one increases abundance and one phosphorylates relationship."""
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
        expected_reified_graph = self.help_make_simple_expected_graph(
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
        """Create a simple reified graph (with 3 nodes - subject, predicate, and object)."""
        expected_reified_graph = nx.DiGraph()
        expected_reified_graph.add_node(edge_num, label=label, causal=causal_tup)
        expected_reified_graph.add_edge(u, edge_num, label=REIF_SUBJECT)
        expected_reified_graph.add_edge(v, edge_num, label=REIF_OBJECT)
        return expected_reified_graph


if __name__ == '__main__':
    unittest.main()
