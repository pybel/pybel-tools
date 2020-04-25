# -*- coding: utf-8 -*-

"""Tests for data integration tools."""

import unittest

from pybel import BELGraph
from pybel.constants import GENE
from pybel.dsl import Gene, Protein, Rna
from pybel_tools.integration import overlay_type_data

HGNC = 'hgnc'


class TestIntegrate(unittest.TestCase):
    """Tests for data integration tools."""

    def test_overlay(self):
        """Test overlaying data in a BEL graph."""
        g = BELGraph()

        g1 = Gene(HGNC, 'a')
        g2 = Gene(HGNC, 'b')
        g3 = Gene(HGNC, 'c')
        g4 = Gene(HGNC, 'd')
        r1 = Rna(HGNC, 'e')
        p1 = Protein(HGNC, 'f')

        g.add_node_from_data(g1)
        g.add_node_from_data(g2)
        g.add_node_from_data(g3)
        g.add_node_from_data(g4)
        g.add_node_from_data(r1)
        g.add_node_from_data(p1)

        self.assertEqual(6, g.number_of_nodes())

        label = 'dgxp'

        overlay_type_data(g, {'a': 1, 'b': 2, 'c': -1}, Gene, HGNC, label=label, impute=0)

        for node in g1, g2, g3, g4:
            self.assertIn(label, g.nodes[node])

        for node in r1, p1:
            self.assertNotIn(label, g.nodes[node])

        self.assertEqual(1, g.nodes[g1][label])
        self.assertEqual(2, g.nodes[g2][label])
        self.assertEqual(-1, g.nodes[g3][label])
        self.assertEqual(0, g.nodes[g4][label])
