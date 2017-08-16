# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel_tools.integration import overlay_type_data

HGNC = 'HGNC'


class TestIntegrate(unittest.TestCase):
    def test_overlay(self):
        g = BELGraph()

        g1 = GENE, HGNC, 'a'
        g2 = GENE, HGNC, 'b'
        g3 = GENE, HGNC, 'c'
        g4 = GENE, HGNC, 'd'
        r1 = RNA, HGNC, 'e'
        p1 = PROTEIN, HGNC, 'f'

        g.add_simple_node(*g1)
        g.add_simple_node(*g2)
        g.add_simple_node(*g3)
        g.add_simple_node(*g4)
        g.add_simple_node(*r1)
        g.add_simple_node(*p1)

        label = 'dgxp'

        overlay_type_data(g, {'a': 1, 'b': 2, 'c': -1}, label, GENE, HGNC, impute=0)

        for node in g1, g2, g3, g4:
            self.assertIn(label, g.node[node])

        for node in r1, p1:
            self.assertNotIn(label, g.node[node])

        self.assertEqual(1, g.node[g1][label])
        self.assertEqual(2, g.node[g2][label])
        self.assertEqual(-1, g.node[g3][label])
        self.assertEqual(0, g.node[g4][label])
