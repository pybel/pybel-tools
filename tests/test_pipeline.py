# -*- coding: utf-8 -*-

import logging
import unittest

from pybel import Pipeline
from pybel.examples import egf_example
from pybel.struct.pipeline import mapped
from pybel.utils import hash_node
from pybel_tools.mutation import build_delete_node_by_hash, build_expand_node_neighborhood_by_hash
from tests.mocks import MockQueryManager

log = logging.getLogger(__name__)
log.setLevel(10)


class TestBoundMutation(unittest.TestCase):
    """Random test for mutation functions"""

    def setUp(self):
        self.graph = egf_example.egf_graph.copy()
        self.original_number_nodes = self.graph.number_of_nodes()
        self.original_number_edges = self.graph.number_of_edges()

        self.manager = MockQueryManager([self.graph])

        build_delete_node_by_hash(self.manager)
        build_expand_node_neighborhood_by_hash(self.manager)

    def tearDown(self):
        """Remove definitions of the built functions"""
        del mapped['expand_node_neighborhood_by_hash']
        del mapped['delete_node_by_hash']

    def check_original_unchanged(self):
        self.assertEqual(self.original_number_nodes, self.graph.number_of_nodes(),
                         msg='original graph nodes should remain unchanged')
        self.assertEqual(self.original_number_edges, self.graph.number_of_edges(),
                         msg='original graph edges should remain unchanged')

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
        self.assertIn('delete_node_by_hash', mapped)
        self.assertIn('expand_node_neighborhood_by_hash', mapped)

    def test_bound_mutation(self):
        """Tests when a node is deleted then re-expanded"""
        pipeline = Pipeline(universe=self.graph)
        pipeline.append('delete_node_by_hash', hash_node(egf_example.nfkb_complex.as_tuple()))
        pipeline.append('expand_node_neighborhood_by_hash', hash_node(egf_example.rela.as_tuple()))

        result = pipeline.run(self.graph, in_place=False)

        self.check_original_unchanged()

        self.assertEqual(self.original_number_nodes, result.number_of_nodes())
        self.assertGreater(self.original_number_edges, result.number_of_edges())
