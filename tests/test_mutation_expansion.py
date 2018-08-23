# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import (
    GENE, RNA,
)
from pybel.dsl import protein, complex_abundance, reaction
from pybel_tools.mutation import enrich_complexes, enrich_reactions

HGNC = 'HGNC'
GOBP = 'GOBP'
CHEBI = 'CHEBI'

g1 = GENE, HGNC, '1'
r1 = RNA, HGNC, '1'
p1_data = protein(HGNC, '1')

g2 = GENE, HGNC, '2'
r2 = RNA, HGNC, '2'
p2_data = protein(HGNC, '2')

complex_1 = complex_abundance(members=[p1_data, p2_data])

reaction_1 = reaction([p1_data], [p2_data])


class TestEnrichGraph(unittest.TestCase):
    def test_enrich_complexes(self):
        graph = BELGraph()

        p2 = graph.add_node_from_data(p2_data)
        c1 = graph.add_node_from_data(complex_1)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        enrich_complexes(graph, graph)

        # The graph should have now one more node (p1) and two new edges connecting the proteins with the complex

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

    def test_enrich_reactions(self):
        graph = BELGraph()

        r1 = graph.add_node_from_data(reaction_1)

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        enrich_reactions(graph, graph)

        # The graph should have now two more nodes (p1, p2) and two new edges connecting the proteins with the reaction

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())
