# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.dsl import abundance, pathology, protein
from pybel_tools.mutation.collapse import collapse_to_protein_interactions

HGNC = 'HGNC'
GOBP = 'GOBP'
CHEBI = 'CHEBI'

g1 = GENE, HGNC, '1'
r1 = RNA, HGNC, '1'
p1_data = protein(HGNC, '1')

g2 = GENE, HGNC, '2'
r2 = RNA, HGNC, '2'
p2_data = protein(HGNC, '2')

g3 = GENE, HGNC, '3'
r3 = RNA, HGNC, '3'
p3 = PROTEIN, HGNC, '3'

g4 = GENE, HGNC, '4'
m4 = MIRNA, HGNC, '4'

a5_data = abundance(CHEBI, '5')
p5_data = pathology(GOBP, '5')


class TestCollapseProteinInteractions(unittest.TestCase):
    def test_protein_interaction_1(self):
        graph = BELGraph()

        p1 = graph.add_node_from_data(p1_data)
        p2 = graph.add_node_from_data(p2_data)
        a5 = graph.add_node_from_data(a5_data)
        p5 = graph.add_node_from_data(p5_data)

        graph.add_edge(p1, p2, **{RELATION: POSITIVE_CORRELATION})
        graph.add_edge(p1, p2, **{RELATION: INCREASES})
        graph.add_edge(a5, p5, **{RELATION: DIRECTLY_INCREASES})
        graph.add_edge(p1, a5, **{RELATION: DECREASES})

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

        graph.add_edge(p1, p2, **{RELATION: POSITIVE_CORRELATION})
        graph.add_edge(p1, p2, **{RELATION: ASSOCIATION})
        graph.add_edge(a5, p5, **{RELATION: DIRECTLY_INCREASES})
        graph.add_edge(p1, a5, **{RELATION: DECREASES})

        collapsed_graph = collapse_to_protein_interactions(graph)

        self.assertEqual(2, collapsed_graph.number_of_nodes())
        self.assertEqual(1, collapsed_graph.number_of_edges())
        self.assertIn(g1, collapsed_graph)
        self.assertIn(g2, collapsed_graph)
