# -*- coding: utf-8 -*-

import logging
import unittest

import pybel_tools.pipeline
import pybel_tools.pipeline
from pybel.examples import egf_example
from pybel.utils import hash_node
from pybel_tools.mutation import build_delete_node_by_hash, build_expand_node_neighborhood_by_hash, infer_central_dogma
from pybel_tools.pipeline import MissingPipelineFunctionError, Pipeline, assert_is_mapped_to_pipeline, mapped
from tests.mocks import MockQueryManager

log = logging.getLogger(__name__)
log.setLevel(10)


class TestEgfExample(unittest.TestCase):
    """Random test for mutation functions"""

    def setUp(self):
        self.graph = egf_example.egf_graph.copy()
        self.original_number_nodes = self.graph.number_of_nodes()
        self.original_number_edges = self.graph.number_of_edges()

    def check_original_unchanged(self):
        self.assertEqual(self.original_number_nodes, self.graph.number_of_nodes(),
                         msg='original graph nodes should remain unchanged')
        self.assertEqual(self.original_number_edges, self.graph.number_of_edges(),
                         msg='original graph edges should remain unchanged')


class TestPipelineFailures(unittest.TestCase):

    def test_assert_failure(self):
        with self.assertRaises(MissingPipelineFunctionError):
            assert_is_mapped_to_pipeline('missing function')

    def test_assert_success(self):
        m = list(mapped)
        self.assertLess(0, len(m))
        m = m[0]
        assert_is_mapped_to_pipeline(m)

    def test_get_function_failure(self):
        pass

    def test_get_function_success(self):
        pass

    def test_fail_add(self):
        pipeline = Pipeline()

        with self.assertRaises(MissingPipelineFunctionError):
            pipeline.append('missing function')

    def test_fail_build(self):
        protocol = [{'function': 'missing function'}]
        with self.assertRaises(MissingPipelineFunctionError):
            Pipeline(protocol=protocol)

    def test_fail_from_json(self):
        protocol = [{'function': 'missing function'}]
        with self.assertRaises(MissingPipelineFunctionError):
            Pipeline.from_json(protocol)


class TestPipeline(TestEgfExample):
    def test_central_dogma_is_registered(self):
        self.assertIn('infer_central_dogma', pybel_tools.pipeline.mapped)

    def test_pipeline_by_string(self):
        pipeline = Pipeline()
        pipeline.append('infer_central_dogma')
        result = pipeline.run(self.graph, in_place=False)

        self.assertEqual(32, result.number_of_nodes())

        for node in self.graph:
            self.assertIn(node, result)

        self.check_original_unchanged()

    def test_pipeline_by_function(self):
        pipeline = Pipeline()
        pipeline.append(infer_central_dogma)
        result = pipeline.run(self.graph, in_place=False)

        self.assertEqual(32, result.number_of_nodes())

        for node in self.graph:
            self.assertIn(node, result)

        self.check_original_unchanged()


class TestBoundMutation(TestEgfExample):
    """Random test for mutation functions"""

    def setUp(self):
        super(TestBoundMutation, self).setUp()

        self.manager = MockQueryManager([self.graph])

        build_delete_node_by_hash(self.manager)
        build_expand_node_neighborhood_by_hash(self.manager)

    def test_mock_contents(self):
        nfkb_complex_tuple = egf_example.nfkb_complex.as_tuple()
        rela_tuple = egf_example.rela.as_tuple()
        self.assertIn(nfkb_complex_tuple, self.manager.graphs[0], msg='Graph missing NFKB complex')
        self.assertIn(rela_tuple, self.manager.graphs[0], msg='Graph missing RELA')

        self.assertIn(nfkb_complex_tuple, self.manager.hash_to_tuple.values(), msg='NFKB is unindexed')
        self.assertIn(rela_tuple, self.manager.hash_to_tuple.values(), msg='RELA is unindexed')

        self.assertIn(hash_node(nfkb_complex_tuple), self.manager.hash_to_tuple.keys(), msg='NFKB is unindexed')
        self.assertIn(hash_node(rela_tuple), self.manager.hash_to_tuple.keys(), msg='RELA is unindexed')

    def test_functions_registered(self):
        self.assertIn('delete_node_by_hash', pybel_tools.pipeline.mapped)
        self.assertIn('expand_node_neighborhood_by_hash', pybel_tools.pipeline.mapped)

    def test_bound_mutation(self):
        """Tests when a node is deleted then re-expanded"""
        pipeline = Pipeline(universe=self.graph)
        pipeline.append('delete_node_by_hash', hash_node(egf_example.nfkb_complex.as_tuple()))
        pipeline.append('expand_node_neighborhood_by_hash', hash_node(egf_example.rela.as_tuple()))

        result = pipeline.run(self.graph, in_place=False)

        self.check_original_unchanged()

        self.assertEqual(self.original_number_nodes, result.number_of_nodes())
        self.assertGreater(self.original_number_edges, result.number_of_edges())
