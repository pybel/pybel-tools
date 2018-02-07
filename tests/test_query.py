# -*- coding: utf-8 -*-

"""
Reference for testing Flask

- Flask Documentation http://flask.pocoo.org/docs/0.12/testing/
- Flask Cookbook: http://flask.pocoo.org/docs/0.12/tutorial/testing/
"""

import logging

from pybel_tools.mutation import collapse_by_central_dogma_to_genes, infer_central_dogma
from pybel_tools.pipeline import Pipeline
from pybel_tools.query import Query
from pybel_tools.selection import get_subgraph_by_annotation_value
from tests.constants import ExampleNetworkMixin, protein_a_tuple, protein_e_tuple, rna_d_tuple
from tests.mocks import MockQueryManager

log = logging.getLogger(__name__)
log.setLevel(10)


class QueryTest(ExampleNetworkMixin):
    def setUp(self):
        super(QueryTest, self).setUp()
        self.manager = MockQueryManager()

    def test_pipeline(self):
        test_network = self.network1  # Defined in test.constants.TestNetworks

        infer_central_dogma(test_network)

        self.assertEqual(9, test_network.number_of_nodes())  # 4 nodes already there +  2*2 proteins + 1 (rna)
        self.assertEqual(8, test_network.number_of_edges())  # 3 already there + 2*2 proteins + 1 (rna)

        network = self.manager.insert_graph(test_network)

        pipeline = Pipeline()
        pipeline.append(collapse_by_central_dogma_to_genes)

        query = Query(
            network_ids=[network.id],
            pipeline=pipeline
        )
        result_graph = query.run(self.manager)

        self.assertEqual(4, result_graph.number_of_nodes())  # same number of nodes than there were
        self.assertEqual(3, result_graph.number_of_edges())  # same number of edges than there were

    def test_pipeline_2(self):
        test_network = self.network1  # Defined in test.constants.TestNetworks

        infer_central_dogma(test_network)

        network = self.manager.insert_graph(test_network)
        network_id = network.id

        pipeline = Pipeline()
        pipeline.append(get_subgraph_by_annotation_value, 'Annotation', 'foo')

        query = Query(network_ids=[network_id])
        query.append_seeding_neighbors([
            protein_a_tuple
        ])
        query.pipeline = pipeline

        result_graph = query.run(self.manager)

        self.assertEqual(3, result_graph.number_of_nodes())  # only expanded to node protein_a and gene_c
        self.assertEqual(2, result_graph.number_of_edges())  # three nodes with two relationships

    def test_query_multiple_networks(self):
        test_network_1 = self.manager.insert_graph(self.network1)
        test_network_2 = self.manager.insert_graph(self.network2)

        query = Query()
        query.append_network(test_network_1.id)
        query.append_network(test_network_2.id)
        query.append_seeding_neighbors([
            rna_d_tuple,
            protein_e_tuple
        ])
        query.append_pipeline(get_subgraph_by_annotation_value, 'Annotation', 'foo')
        query.append_pipeline(collapse_by_central_dogma_to_genes)

        result_graph = query.run(self.manager)

        self.assertEqual(4, result_graph.number_of_nodes())
        # TODO: discuss this with Charlie. It would be cool to infer the edge between b and a
        self.assertEqual(2, result_graph.number_of_edges())

    def test_query_multiple_networks_with_api(self):
        test_network_1 = self.manager.insert_graph(self.network1)
        test_network_2 = self.manager.insert_graph(self.network2)

        pipeline = Pipeline()
        pipeline.append(get_subgraph_by_annotation_value, 'Annotation', 'foo')
        pipeline.append(collapse_by_central_dogma_to_genes)

        query = Query(
            network_ids=[test_network_1.id, test_network_2.id],
            pipeline=pipeline
        )
        query.append_seeding_neighbors([rna_d_tuple, protein_e_tuple])

        result_graph = query.run(self.manager)

        self.assertEqual(4, result_graph.number_of_nodes())
        # TODO: discuss this with Charlie. It would be cool to infer the edge between b and a
        self.assertEqual(2, result_graph.number_of_edges())
