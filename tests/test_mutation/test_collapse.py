# -*- coding: utf-8 -*-

"""Tests for collapse functions."""

import unittest

from pybel import BELGraph
from pybel.constants import (
    ASSOCIATION, DECREASES, DIRECTLY_INCREASES, INCREASES, POSITIVE_CORRELATION,
)
from pybel.dsl import abundance, gene, mirna, pathology, protein, rna
from pybel.testing.utils import n
from pybel_tools.mutation.collapse import collapse_to_protein_interactions

HGNC = 'HGNC'
GOBP = 'GOBP'
CHEBI = 'CHEBI'

g1 = gene(HGNC, '1')
r1 = rna(HGNC, '1')
p1_data = protein(HGNC, '1')

g2 = gene(HGNC, '2')
r2 = rna(HGNC, '2')
p2_data = protein(HGNC, '2')

g3 = gene(HGNC, '3')
r3 = rna(HGNC, '3')
p3 = protein(HGNC, '3')

g4 = gene(HGNC, '4')
m4 = mirna(HGNC, '4')

a5_data = abundance(CHEBI, '5')
p5_data = pathology(GOBP, '5')


class TestCollapseProteinInteractions(unittest.TestCase):
    def test_protein_interaction_1(self):
        graph = BELGraph()

        p1 = graph.add_node_from_data(p1_data)
        p2 = graph.add_node_from_data(p2_data)
        a5 = graph.add_node_from_data(a5_data)
        p5 = graph.add_node_from_data(p5_data)

        graph.add_qualified_edge(p1, p2, POSITIVE_CORRELATION, n(), n())

        graph.add_qualified_edge(p1, p2, INCREASES, n(), n())
        graph.add_qualified_edge(a5, p5, DIRECTLY_INCREASES, n(), n())
        graph.add_qualified_edge(p1, a5, DECREASES, n(), n())

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

        graph.add_qualified_edge(p1, p2, POSITIVE_CORRELATION, n(), n())
        graph.add_qualified_edge(p1, p2, ASSOCIATION, n(), n())
        graph.add_qualified_edge(a5, p5, DIRECTLY_INCREASES, n(), n())
        graph.add_qualified_edge(p1, a5, DECREASES, n(), n())

        collapsed_graph = collapse_to_protein_interactions(graph)

        self.assertEqual(2, collapsed_graph.number_of_nodes())
        self.assertEqual(1, collapsed_graph.number_of_edges())
        self.assertIn(g1, collapsed_graph)
        self.assertIn(g2, collapsed_graph)
