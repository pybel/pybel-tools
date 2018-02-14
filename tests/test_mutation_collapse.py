# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.constants import unqualified_edge_code
from pybel_tools.mutation import collapse_by_central_dogma, collapse_nodes
from pybel_tools.mutation.collapse import collapse_to_protein_interactions

HGNC = 'HGNC'
GOBP = 'GOBP'
CHEBI = 'CHEBI'

g1 = GENE, HGNC, '1'
r1 = RNA, HGNC, '1'
p1 = PROTEIN, HGNC, '1'

g2 = GENE, HGNC, '2'
r2 = RNA, HGNC, '2'
p2 = PROTEIN, HGNC, '2'

g3 = GENE, HGNC, '3'
r3 = RNA, HGNC, '3'
p3 = PROTEIN, HGNC, '3'

g4 = GENE, HGNC, '4'
m4 = MIRNA, HGNC, '4'

a5 = ABUNDANCE, CHEBI, '5'
p5 = PATHOLOGY, GOBP, '5'


class TestCollapseProteinInteractions(unittest.TestCase):
    def test_protein_interaction_1(self):
        graph = BELGraph()

        graph.add_simple_node(*p1)
        graph.add_simple_node(*p2)
        graph.add_simple_node(*a5)
        graph.add_simple_node(*p5)

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

        graph.add_simple_node(*p1)
        graph.add_simple_node(*p2)
        graph.add_simple_node(*a5)
        graph.add_simple_node(*p5)

        graph.add_edge(p1, p2, **{RELATION: POSITIVE_CORRELATION})
        graph.add_edge(p1, p2, **{RELATION: ASSOCIATION})
        graph.add_edge(a5, p5, **{RELATION: DIRECTLY_INCREASES})
        graph.add_edge(p1, a5, **{RELATION: DECREASES})

        collapsed_graph = collapse_to_protein_interactions(graph)

        self.assertEqual(2, collapsed_graph.number_of_nodes())
        self.assertEqual(1, collapsed_graph.number_of_edges())
        self.assertIn(g1, collapsed_graph)
        self.assertIn(g2, collapsed_graph)


class TestCollapseDownstream(unittest.TestCase):
    def test_collapse_1(self):
        graph = BELGraph()

        graph.add_simple_node(*p1)
        graph.add_simple_node(*p2)
        graph.add_simple_node(*p3)

        graph.add_edge(p1, p3, **{RELATION: INCREASES})
        graph.add_edge(p2, p3, **{RELATION: DIRECTLY_INCREASES})

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        d = {
            p1: {p2}
        }

        collapse_nodes(graph, d)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges(), msg=graph.edges(data=True, keys=True))

    def test_collapse_dogma_1(self):
        graph = BELGraph()

        graph.add_simple_node(*p1)
        graph.add_simple_node(*r1)

        graph.add_edge(r1, p1, key=unqualified_edge_code[TRANSLATED_TO], **{RELATION: TRANSLATED_TO})

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        collapse_by_central_dogma(graph)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

    def test_collapse_dogma_2(self):
        graph = BELGraph()

        graph.add_simple_node(*p1)
        graph.add_simple_node(*r1)
        graph.add_simple_node(*g1)

        graph.add_edge(r1, p1, key=unqualified_edge_code[TRANSLATED_TO], **{RELATION: TRANSLATED_TO})
        graph.add_edge(g1, r1, key=unqualified_edge_code[TRANSCRIBED_TO], **{RELATION: TRANSCRIBED_TO})

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        collapse_by_central_dogma(graph)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

    def test_collapse_dogma_3(self):
        graph = BELGraph()

        graph.add_simple_node(*r1)
        graph.add_simple_node(*g1)

        graph.add_edge(g1, r1, key=unqualified_edge_code[TRANSCRIBED_TO], **{RELATION: TRANSCRIBED_TO})

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        collapse_by_central_dogma(graph)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())
