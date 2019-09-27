# -*- coding: utf-8 -*-

import unittest

import pybel
from pybel.dsl import bioprocess, protein
from pybel.testing.utils import n
from pybel_tools.generation import generate_bioprocess_mechanisms


class TestGenerate(unittest.TestCase):
    def test_simple(self):
        graph = pybel.BELGraph()

        key = 'DGXP'

        a = protein('HGNC', 'A')
        b = protein('HGNC', 'B')
        c = protein('HGNC', 'c')
        d = bioprocess('GOBP', 'D')

        graph.add_node_from_data(a)
        graph.add_node_from_data(b)
        graph.add_node_from_data(c)
        graph.add_node_from_data(d)

        graph.nodes[a][key] = 2
        graph.nodes[b][key] = -1
        graph.nodes[c][key] = 1

        graph.add_increases(a, b, citation=n(), evidence=n())
        graph.add_decreases(b, d, citation=n(), evidence=n())
        graph.add_increases(a, c, citation=n(), evidence=n())
        graph.add_increases(c, d, citation=n(), evidence=n())

        candidate_mechanisms = generate_bioprocess_mechanisms(graph, key)

        self.assertEqual(1, len(candidate_mechanisms))
        self.assertIn(d, candidate_mechanisms)

        # score = heat.workflow_average(graph, d, key, runs=5)
        # self.assertEqual(3, score)


if __name__ == '__main__':
    unittest.main()
