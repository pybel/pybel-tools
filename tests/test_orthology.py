# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel_tools.orthology import HGNC, MGI
from pybel_tools.orthology import integrate_orthologies_from_rgd, collapse_orthologies
from tests.constants import rgd_orthologs_path


class TestOrthology(unittest.TestCase):
    def test_collapse_by_orthology(self):
        graph = BELGraph()

        g1 = GENE, HGNC, 'A1BG'
        g2 = GENE, HGNC, 'b'
        g3 = GENE, MGI, 'A1bg'
        g4 = GENE, MGI, 'c'

        graph.add_simple_node(*g1)
        graph.add_simple_node(*g2)
        graph.add_simple_node(*g3)
        graph.add_simple_node(*g4)

        self.assertEqual(4, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        integrate_orthologies_from_rgd(graph, rgd_orthologs_path)

        self.assertEqual(4, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        self.assertTrue(graph.has_edge(g1, g3))
        self.assertFalse(graph.has_edge(g3, g1))

        collapse_orthologies(graph)

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())
