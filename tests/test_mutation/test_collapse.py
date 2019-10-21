# -*- coding: utf-8 -*-

"""Tests for collapse functions."""

import unittest

from pybel import BELGraph
from pybel.constants import (
    ASSOCIATION, DECREASES, DIRECTLY_INCREASES, INCREASES, POSITIVE_CORRELATION,
)
from pybel.dsl import Gene, Protein, abundance, mirna, pathology, rna
from pybel.testing.utils import n

from pybel_tools.mutation.collapse import collapse_to_protein_interactions

HGNC = 'HGNC'
GO = 'GO'
CHEBI = 'CHEBI'

g1 = Gene(HGNC, '1')
r1 = rna(HGNC, '1')
p1_data = Protein(HGNC, '1')

g2 = Gene(HGNC, '2')
r2 = rna(HGNC, '2')
p2_data = Protein(HGNC, '2')

g3 = Gene(HGNC, '3')
r3 = rna(HGNC, '3')
p3 = Protein(HGNC, '3')

g4 = Gene(HGNC, '4')
m4 = mirna(HGNC, '4')

a5_data = abundance(CHEBI, '5')
p5_data = pathology(GO, '5')


class TestCollapseProteinInteractions(unittest.TestCase):
    def test_protein_interaction_1(self):
        graph = BELGraph()

        p1 = graph.add_node_from_data(p1_data)
        p2 = graph.add_node_from_data(p2_data)
        a5 = graph.add_node_from_data(a5_data)
        p5 = graph.add_node_from_data(p5_data)

        graph.add_qualified_edge(p1, p2, relation=POSITIVE_CORRELATION, citation=n(), evidence=n())

        graph.add_qualified_edge(p1, p2, relation=INCREASES, citation=n(), evidence=n())
        graph.add_qualified_edge(a5, p5, relation=DIRECTLY_INCREASES, citation=n(), evidence=n())
        graph.add_qualified_edge(p1, a5, relation=DECREASES, citation=n(), evidence=n())

        collapsed_graph = collapse_to_protein_interactions(graph)

        self.assertEqual(2, collapsed_graph.number_of_nodes())
        self.assertEqual(2, collapsed_graph.number_of_edges())
        self.assertIn(g1, collapsed_graph)
        self.assertIn(g2, collapsed_graph)

    def test_protein_interaction_2(self):
        graph = BELGraph()

        p1 = graph.add_node_from_data(p1_data)
        p2 = graph.add_node_from_data(p2_data)
        a5 = graph.add_node_from_data(a5_data)
        p5 = graph.add_node_from_data(p5_data)

        graph.add_qualified_edge(p1, p2, relation=POSITIVE_CORRELATION, citation=n(), evidence=n())
        graph.add_qualified_edge(p1, p2, relation=ASSOCIATION, citation=n(), evidence=n())
        graph.add_qualified_edge(a5, p5, relation=DIRECTLY_INCREASES, citation=n(), evidence=n())
        graph.add_qualified_edge(p1, a5, relation=DECREASES, citation=n(), evidence=n())

        collapsed_graph = collapse_to_protein_interactions(graph)

        self.assertEqual(2, collapsed_graph.number_of_nodes())
        self.assertEqual(1, collapsed_graph.number_of_edges())
        self.assertIn(g1, collapsed_graph)
        self.assertIn(g2, collapsed_graph)
