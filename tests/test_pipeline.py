# -*- coding: utf-8 -*-

import logging
import unittest

import pybel_tools.pipeline
import pybel_tools.pipeline
from pybel.examples import egf_example
from pybel_tools.mutation import build_delete_node_by_hash, build_expand_node_neighborhood_by_hash, infer_central_dogma
from pybel_tools.pipeline import Pipeline
from tests.mocks import MockQueryManager

log = logging.getLogger(__name__)
log.setLevel(10)


class TestPipeline(unittest.TestCase):
    def setUp(self):
        self.graph = egf_example.egf_graph
        self.original_number_nodes = self.graph.number_of_nodes()

        self.original_number_edges = self.graph.number_of_edges()

    def check_original_unchanged(self):
        self.assertEqual(self.original_number_nodes, self.graph.number_of_nodes(),
                         msg='original graph nodes should remain unchanged')
        self.assertEqual(self.original_number_edges, self.graph.number_of_edges(),
                         msg='original graph edges should remain unchanged')

    def test_infer_central_dogma_lookup(self):
        func = 'infer_central_dogma'
        self.assertIn(func, pybel_tools.pipeline.mapped)

        pipeline = Pipeline()
        pipeline.append(func)
        result = pipeline.run(self.graph)
        self.check_original_unchanged()
        self.assertEqual(12 + 10 * 2, result.number_of_nodes())

    def test_infer_central_dogma(self):
        pipeline = Pipeline()
        pipeline.append(infer_central_dogma)
        result = pipeline.run(self.graph)
        self.check_original_unchanged()
        self.assertEqual(12 + 10 * 2, result.number_of_nodes())


class TestBoundMutation(unittest.TestCase):
    """Random test for mutation functions"""

    def test_bound_mutation(self):
        """Tests when a node is deleted then re-expanded"""
        graph = egf_example.egf_graph

        original_number_nodes = graph.number_of_nodes()
        original_number_edges = graph.number_of_edges()

        manager = MockQueryManager([graph])

        self.assertIn(graph.hash_node(egf_example.nfkb_complex), manager.hash_to_tuple)
        self.assertIn(graph.hash_node(egf_example.rela), manager.hash_to_tuple)

        build_delete_node_by_hash(manager)
        self.assertIn('delete_node_by_hash', pybel_tools.pipeline.mapped)

        build_expand_node_neighborhood_by_hash(manager)
        self.assertIn('expand_node_neighborhood_by_hash', pybel_tools.pipeline.mapped)

        pipeline = Pipeline(universe=graph)
        pipeline.append('delete_node_by_hash', graph.hash_node(egf_example.nfkb_complex))
        pipeline.append('expand_node_neighborhood_by_hash', graph.hash_node(egf_example.rela))

        result = pipeline.run(graph)

        self.assertEqual(original_number_nodes, graph.number_of_nodes(),
                         msg='original graph nodes should remain unchanged')
        self.assertEqual(original_number_edges, graph.number_of_edges(),
                         msg='original graph edges should remain unchanged')

        self.assertEqual(original_number_nodes, result.number_of_nodes())
        self.assertLess(graph.number_of_edges(), original_number_edges)
