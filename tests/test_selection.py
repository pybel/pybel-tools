# -*- coding: utf-8 -*-

import itertools as itt
import random
import unittest
from uuid import uuid4

import numpy as np

from pybel import BELGraph
from pybel.constants import INCREASES
from pybel.dsl import protein
from pybel_tools.selection import get_random_subgraph


class TestRandomSample(unittest.TestCase):
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
