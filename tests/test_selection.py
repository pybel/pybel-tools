# -*- coding: utf-8 -*-

import itertools as itt
import random
import unittest
from collections import Counter
from uuid import uuid4

import numpy as np

from pybel import BELGraph
from pybel.constants import INCREASES
from pybel.dsl import protein
from pybel_tools.selection import get_random_subgraph
from pybel_tools.selection.induce_subgraph import randomly_select_node


class TestRandomSample(unittest.TestCase):
    def test_randomly_select_node_1(self):
        g = BELGraph()
        g.add_edge(1, 2)
        g.add_edge(2, 3)
        g.add_edge(4, 3)

        no_grow = {3}
        random_state = np.random.RandomState(seed=127)

        nodes = []
        trials = 1000
        for _ in range(trials):
            node = randomly_select_node(g, no_grow, random_state)
            nodes.append(node)

        node_counter = Counter(nodes)

        self.assertAlmostEqual(.1, node_counter[1] / trials, places=1)
        self.assertAlmostEqual(.2, node_counter[2] / trials, places=1)
        self.assertAlmostEqual(.2, node_counter[3] / trials, places=1)
        self.assertAlmostEqual(.3, node_counter[4] / trials, places=1)

    def test_okay(self):
        np.random.seed(127)

        g = BELGraph()

        n_edges = 500

        nodes = [
            protein(namespace='NS', name=str(i))
            for i in range(1, 50)
        ]

        edges = list(itt.combinations(nodes, r=2))
        random.shuffle(edges)

        for u, v in edges[:n_edges]:
            g.add_qualified_edge(
                u, v,
                relation=INCREASES,
                citation=str(uuid4()),
                evidence=str(uuid4()),
            )

        self.assertEqual(n_edges, g.number_of_edges())

        sg = get_random_subgraph(g, number_edges=250, number_seed_nodes=5, seed=127)
        self.assertEqual(250, sg.number_of_edges())

    def test_too_small(self):
        np.random.seed(127)

        g = BELGraph()

        n_edges = 25

        nodes = [
            protein(namespace='NS', name=str(i))
            for i in range(1, 11)
        ]

        edges = list(itt.combinations(nodes, r=2))
        random.shuffle(edges)

        for u, v in edges[:n_edges]:
            print(u, v)
            g.add_qualified_edge(
                u, v,
                relation=INCREASES,
                citation=str(uuid4()),
                evidence=str(uuid4()),
            )

        self.assertEqual(n_edges, g.number_of_edges())

        sg = get_random_subgraph(g, number_edges=250, number_seed_nodes=5, seed=127)

        self.assertEqual(g.number_of_edges(), sg.number_of_edges(),
                         msg='since graph is too small, the subgraph should contain the whole thing')
