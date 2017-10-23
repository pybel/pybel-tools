# -*- coding: utf-8 -*-

"""
Reference for testing Flask

- Flask Documentation http://flask.pocoo.org/docs/0.12/testing/
- Flask Cookbook: http://flask.pocoo.org/docs/0.12/tutorial/testing/
"""

import logging

import networkx as nx

from pybel.constants import *
from pybel.utils import hash_node
from pybel_tools.mutation import *
from pybel_tools.pipeline import Pipeline
from pybel_tools.query import Query
from pybel_tools.selection import *
from tests.constants import ExampleNetworkMixin, ManagerMixin

HGNC = 'HGNC'

log = logging.getLogger(__name__)
log.setLevel(10)


def relabel_nodes_to_hashes(graph, copy=False):
    """Relabels nodes to hashes in the graph

    :param pybel.BELGraph graph:
    :param bool copy: Copy the graph?
    :rtype: pybel.BELGraph
    """
    if 'PYBEL_RELABELED' in graph.graph:
        log.warning('%s has already been relabeled', graph.name)
        return graph

    mapping = {}
    for node in graph:
        mapping[node] = hash_node(node)

    nx.relabel.relabel_nodes(graph, mapping, copy=copy)

    graph.graph['PYBEL_RELABELED'] = True

    return graph


class QueryTest(ExampleNetworkMixin, ManagerMixin):
    def test_pipeline(self):
        test_network = self.network1  # Defined in test.constants.TestNetworks

        infer_central_dogma(test_network)

        self.assertEqual(9, test_network.number_of_nodes())  # 4 nodes already there +  2*2 proteins + 1 (rna)
        self.assertEqual(8, test_network.number_of_edges())  # 3 already there + 2*2 proteins + 1 (rna)

        network = self.manager.insert_graph(test_network, store_parts=False)

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

        network = self.manager.insert_graph(test_network, store_parts=False)
        network_id = network.id

        pipeline = Pipeline()
        pipeline.append(get_subgraph_by_annotation_value, 'Annotation', 'foo')

        query = Query(
            network_ids=network_id,
            seed_list=[{'type': 'neighbors', 'data': [(PROTEIN, HGNC, 'a')]}],
            pipeline=pipeline
        )

        result_graph = query.run(self.manager)

        self.assertEqual(3, result_graph.number_of_nodes())  # only expanded to node protein_a and gene_c
        self.assertEqual(2, result_graph.number_of_edges())  # three nodes with two relationships

    def test_query_multiple_networks(self):
        test_network_1 = self.manager.insert_graph(self.network1)
        test_network_2 = self.manager.insert_graph(self.network2)

        query = Query(network_ids=[test_network_1.id, test_network_2.id])
        query.add_seed_neighbors([
            (RNA, HGNC, 'd'),
            (PROTEIN, HGNC, 'e')
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

        node_1 = RNA, HGNC, 'd'
        node_2 = PROTEIN, HGNC, 'e'

        # node_1_id = api.get_node_id(node_1)
        # node_2_id = api.get_node_id(node_2)

        query = Query(
            network_ids=[test_network_1.id, test_network_2.id],
            pipeline=pipeline
        )
        query.add_seed_neighbors([node_1, node_2])

        result_graph = query.run(self.manager)
        result_graph = relabel_nodes_to_hashes(result_graph)

        self.assertEqual(4, result_graph.number_of_nodes())
        # TODO: discuss this with Charlie. It would be cool to infer the edge between b and a
        self.assertEqual(2, result_graph.number_of_edges())

    def test_get_network_by_name(self):
        test_network_1 = self.network1  # Defined in test.constants.TestNetworks
        test_network_2 = self.network2  # Defined in test.constants.TestNetworks

        test_network_2 = self.manager.insert_graph(test_network_2, store_parts=False)
        test_network_1 = self.manager.insert_graph(test_network_1, store_parts=False)  # Latest updated network

        network = self.manager.get_most_recent_network_by_name('network_test')

        self.assertEqual('1.1.0', network.version)
        self.assertEqual(test_network_1.version, network.version)
