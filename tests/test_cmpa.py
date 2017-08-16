# -*- coding: utf-8 -*-

import unittest

import pybel
from pybel.constants import PROTEIN, BIOPROCESS, RELATION, INCREASES, DECREASES
from pybel_tools.analysis import cmpa


class TestCmpa(unittest.TestCase):
    def test_simple(self):
        graph = pybel.BELGraph()

        key = 'DGXP'

        a = PROTEIN, 'HGNC', 'A'
        b = PROTEIN, 'HGNC', 'B'
        c = PROTEIN, 'HGNC', 'c'
        d = BIOPROCESS, 'GOBP', 'D'

        graph.add_simple_node(*a)
        graph.add_simple_node(*b)
        graph.add_simple_node(*c)
        graph.add_simple_node(*d)

        graph.node[a][key] = 2
        graph.node[b][key] = -1
        graph.node[c][key] = 1

        graph.add_edge(a, b, attr_dict={RELATION: INCREASES})
        graph.add_edge(b, d, attr_dict={RELATION: DECREASES})
        graph.add_edge(a, c, attr_dict={RELATION: INCREASES})
        graph.add_edge(c, d, attr_dict={RELATION: INCREASES})

        candidate_mechanisms = cmpa.generate_bioprocess_mechanisms(graph, key)

        self.assertEqual(1, len(candidate_mechanisms))
        self.assertIn(d, candidate_mechanisms)

        # score = cmpa.workflow_average(graph, d, key, runs=5)
        # self.assertEqual(3, score)


if __name__ == '__main__':
    unittest.main()
