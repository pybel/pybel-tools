# -*- coding: utf-8 -*-

"""Tests for node utilities."""

import unittest

from pybel import BELGraph
from pybel.dsl import Protein
from pybel.testing.utils import n
from pybel_tools.mutation.collapse import collapse_nodes_with_same_names


class TestMerge(unittest.TestCase):
    def test_merge_nodes_by_name(self):
        graph = BELGraph()
        priority = ['ncbigene', 'hgnc']
        a_hgnc, a_entrez = [Protein(namespace=namespace, name='a') for namespace in ('hgnc', 'ncbigene')]
        b = Protein('ncbigene', 'b')
        graph.add_increases(a_hgnc, b, citation=n(), evidence=n())
        graph.add_increases(a_entrez, b, citation=n(), evidence=n())

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        collapse_nodes_with_same_names(graph, priority=priority, use_tqdm=False)
        self.assertEqual(2, graph.number_of_nodes(), msg='Wrong number remaining nodes: {}'.format(graph.nodes()))
        self.assertEqual(2, graph.number_of_edges(), msg=f'Wrong number remaining edges: {graph.edges()}')
        self.assertIn(a_entrez, graph, msg=f'Nodes: {list(graph)}')
        self.assertIn(b, graph)
        self.assertNotIn(a_hgnc, graph)
