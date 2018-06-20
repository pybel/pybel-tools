# -*- coding: utf-8 -*-

import unittest

import pybel
from pybel.constants import DECREASES, INCREASES
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

        graph.node[a.as_tuple()][key] = 2
        graph.node[b.as_tuple()][key] = -1
        graph.node[c.as_tuple()][key] = 1

        graph.add_qualified_edge(a, b, INCREASES, n(), n())
        graph.add_qualified_edge(b, d, DECREASES, n(), n())
        graph.add_qualified_edge(a, c, INCREASES, n(), n())
        graph.add_qualified_edge(c, d, INCREASES, n(), n())

        candidate_mechanisms = generate_bioprocess_mechanisms(graph, key)

        self.assertEqual(1, len(candidate_mechanisms))
        self.assertIn(d.as_tuple(), candidate_mechanisms)

        # score = heat.workflow_average(graph, d, key, runs=5)
        # self.assertEqual(3, score)


if __name__ == '__main__':
    unittest.main()
