# -*- coding: utf-8 -*-

import itertools as itt
import random
import unittest
from collections import Counter
from uuid import uuid4

import numpy as np

from pybel import BELGraph
from pybel.constants import INCREASES, PROTEIN
from pybel.dsl import protein
from pybel_tools.selection import get_random_subgraph
from pybel_tools.selection.random_subgraph import randomly_select_node


def n():
    """Generates a PyBEL node tuple

    :rtype: tuple
    """
    return PROTEIN, 'TEST', str(uuid4())


class TestRandomSelectNode(unittest.TestCase):
    """Test random node selection"""

    def setUp(self):
        self.random_state = np.random.RandomState(seed=127)
        self.trials = 30000

    def test_randomly_select_node_1(self):
        """Tests that randomly selecting nodes works"""
        a, b, c, d = (n() for _ in range(4))

        g = BELGraph()
        g.add_edge(a, b)
        g.add_edge(b, c)
        g.add_edge(b, d)

        self.assertEqual(1, g.degree(a))
        self.assertEqual(3, g.degree(b))
        self.assertEqual(1, g.degree(c))
        self.assertEqual(1, g.degree(d))

        no_grow = set()

        node_counter = Counter(
            randomly_select_node(g, no_grow, self.random_state)
            for _ in range(self.trials)
        )

        self.assertIn(a, node_counter)
        self.assertAlmostEqual((1 / 6), node_counter[a] / self.trials, places=2)

        self.assertIn(b, node_counter)
        self.assertAlmostEqual((3 / 6), node_counter[b] / self.trials, places=2)

        self.assertIn(c, node_counter)
        self.assertAlmostEqual((1 / 6), node_counter[c] / self.trials, places=2)

        self.assertIn(d, node_counter)
        self.assertAlmostEqual((1 / 6), node_counter[d] / self.trials, places=2)

    def test_randomly_select_node_2(self):
        """Tests that randomly selecting nodes works, but disallow C"""
        a, b, c, d = (n() for _ in range(4))

        g = BELGraph()
        g.add_edge(a, b)
        g.add_edge(b, c)
        g.add_edge(b, d)

        self.assertEqual(1, g.degree(a))
        self.assertEqual(3, g.degree(b))
        self.assertEqual(1, g.degree(c))
        self.assertEqual(1, g.degree(d))

        no_grow = {c}

        node_counter = Counter(
            randomly_select_node(g, no_grow, self.random_state)
            for _ in range(self.trials)
        )

        self.assertIn(a, node_counter)
        self.assertAlmostEqual((1 / 5), node_counter[a] / self.trials, places=2)

        self.assertIn(b, node_counter)
        self.assertAlmostEqual((3 / 5), node_counter[b] / self.trials, places=2)

        self.assertNotIn(c, node_counter)

        self.assertIn(d, node_counter)
        self.assertAlmostEqual((1 / 5), node_counter[d] / self.trials, places=2)


def make_nodes(n):
    """Returns a list of PyBEL node data dictionaries

    :param int n: number nodes
    :rtype: list[protein]
    """
    return [
        protein(namespace='NS', name=str(i))
        for i in range(1, n)
    ]


class TestRandomSample(unittest.TestCase):
    def setUp(self):
        self.random_state = np.random.seed(127)
        random.seed(127)

    def test_okay(self):
        graph = BELGraph()
        nodes = make_nodes(50)

        edges = list(itt.combinations(nodes, r=2))
        random.shuffle(edges)

        n_edges = 500

        for u, v in edges[:n_edges]:
            graph.add_qualified_edge(
                u, v,
                relation=INCREASES,
                citation=str(uuid4()),
                evidence=str(uuid4()),
            )

        self.assertEqual(n_edges, graph.number_of_edges())

        sg = get_random_subgraph(graph, number_edges=250, number_seed_edges=5, seed=127)
        self.assertEqual(250, sg.number_of_edges())

    def test_too_small(self):
        graph = BELGraph()
        nodes = make_nodes(11)

        edges = list(itt.combinations(nodes, r=2))
        random.shuffle(edges)

        n_edges = 25

        for u, v in edges[:n_edges]:
            graph.add_qualified_edge(
                u, v,
                relation=INCREASES,
                citation=str(uuid4()),
                evidence=str(uuid4()),
            )

        self.assertEqual(n_edges, graph.number_of_edges())

        sg = get_random_subgraph(graph, number_edges=250, number_seed_edges=5, seed=127)

        self.assertEqual(graph.number_of_edges(), sg.number_of_edges(),
                         msg='since graph is too small, the subgraph should contain the whole thing')
