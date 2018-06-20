# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.dsl import gene, protein, rna
from pybel_tools.integration import overlay_type_data

HGNC = 'HGNC'


class TestIntegrate(unittest.TestCase):
    def test_overlay(self):
        g = BELGraph()

        g1 = gene(HGNC, 'a')
        g2 = gene(HGNC, 'b')
        g3 = gene(HGNC, 'c')
        g4 = gene(HGNC, 'd')
        r1 = rna(HGNC, 'e')
        p1 = protein(HGNC, 'f')

        g.add_node_from_data(g1)
        g.add_node_from_data(g2)
        g.add_node_from_data(g3)
        g.add_node_from_data(g4)
        g.add_node_from_data(r1)
        g.add_node_from_data(p1)

        self.assertEqual(6, g.number_of_nodes())

        label = 'dgxp'

        overlay_type_data(g, {'a': 1, 'b': 2, 'c': -1}, GENE, HGNC, label=label, impute=0)

        for node in g1, g2, g3, g4:
            self.assertIn(label, g.node[node.as_tuple()])

        for node in r1, p1:
            self.assertNotIn(label, g.node[node.as_tuple()])

        self.assertEqual(1, g.node[g1.as_tuple()][label])
        self.assertEqual(2, g.node[g2.as_tuple()][label])
        self.assertEqual(-1, g.node[g3.as_tuple()][label])
        self.assertEqual(0, g.node[g4.as_tuple()][label])
