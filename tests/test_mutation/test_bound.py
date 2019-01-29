# -*- coding: utf-8 -*-

"""Tests for bound mutation functions."""

import unittest

from pybel import Pipeline
from pybel.examples.egf_example import egf_graph, nfkb_complex, rela
from pybel.struct.pipeline import mapped
from pybel.testing.mock_manager import MockQueryManager
from pybel_tools.mutation import build_delete_node_by_hash, build_expand_node_neighborhood_by_hash


class TestBoundMutation(unittest.TestCase):
    """Random test for mutation functions."""

    def setUp(self):
        self.graph = egf_graph.copy()
        self.original_number_nodes = self.graph.number_of_nodes()
        self.original_number_edges = self.graph.number_of_edges()

        self.manager = MockQueryManager()
        self.network_id = self.manager.insert_graph(self.graph).id

        build_delete_node_by_hash(self.manager)
        build_expand_node_neighborhood_by_hash(self.manager)

    def tearDown(self):
        """Remove definitions of the built functions"""
        del mapped['expand_node_neighborhood_by_hash']
        del mapped['delete_node_by_hash']

    def test_functions_registered(self):
        self.assertIn('delete_node_by_hash', mapped)
        self.assertIn('expand_node_neighborhood_by_hash', mapped)

    def check_original_unchanged(self):
        self.assertEqual(self.original_number_nodes, self.graph.number_of_nodes(),
                         msg='original graph nodes should remain unchanged')
        self.assertEqual(self.original_number_edges, self.graph.number_of_edges(),
                         msg='original graph edges should remain unchanged')

    def test_mock_contents(self):
        self.assertIn(nfkb_complex, self.graph, msg='Graph missing NFKB complex')
        self.assertIn(rela, self.graph, msg='Graph missing RELA')

        self.assertIn(nfkb_complex.sha512, self.manager.hash_to_node, msg='NFKB is unindexed')
        self.assertIn(rela.sha512, self.manager.hash_to_node, msg='RELA is unindexed')

        self.assertIn(nfkb_complex, self.manager.hash_to_node.values(), msg='NFKB is unindexed')
        self.assertIn(rela, self.manager.hash_to_node.values(), msg='RELA is unindexed')

    def test_bound_mutation(self):
        """Test when a node is deleted then re-expanded."""
        pipeline = Pipeline()
        pipeline.append('delete_node_by_hash', nfkb_complex.sha512)
        pipeline.append('expand_node_neighborhood_by_hash', rela.sha512)
        result = pipeline.run(self.graph)

        self.check_original_unchanged()

        self.assertEqual(self.original_number_nodes, result.number_of_nodes())
        self.assertGreater(self.original_number_edges, result.number_of_edges())
