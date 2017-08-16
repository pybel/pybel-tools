# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel_tools.mutation import expand_node_neighborhood, expand_all_node_neighborhoods


class TestExpandNeighborhood(unittest.TestCase):
    def setUp(self):
        self.universe = BELGraph()

        self.universe.add_edge(1, 2)
        self.universe.add_edge(2, 3)
        self.universe.add_edge(3, 7)
        self.universe.add_edge(1, 4)
        self.universe.add_edge(1, 5)
        self.universe.add_edge(5, 6)
        self.universe.add_edge(8, 2)

        self.graph = BELGraph()
        self.graph.add_edge(1, 2)

    def test_expand_failure(self):
        self.graph.add_node(0)
        with self.assertRaises(Exception):
            expand_node_neighborhood(self.universe, self.graph, 0)

    def test_expand_add(self):
        self.assertNotIn(3, self.graph)

        expand_node_neighborhood(self.universe, self.graph, 3)

        self.assertIn(3, self.graph)
        self.assertIn(7, self.graph)
        self.assertIn(7, self.graph.edge[3])

    def test_expand_successors(self):
        expand_node_neighborhood(self.universe, self.graph, 1)

        self.assertIn(4, self.graph)
        self.assertIn(5, self.graph)
        self.assertIn(5, self.graph.edge[1])

    def test_expand_predecessors(self):
        expand_node_neighborhood(self.universe, self.graph, 2)

        self.assertIn(8, self.graph)
        self.assertIn(2, self.graph.edge[8])

    def test_expand_all_neighborhoods(self):
        expand_all_node_neighborhoods(self.universe, self.graph)
        self.assertIn(3, self.graph)
        self.assertIn(3, self.graph.edge[2])

        self.assertIn(4, self.graph)
        self.assertIn(4, self.graph.edge[1])
        self.assertIn(5, self.graph)
        self.assertIn(5, self.graph.edge[1])

        self.assertIn(8, self.graph)
        self.assertIn(2, self.graph.edge[8])
