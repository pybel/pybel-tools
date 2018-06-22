# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.dsl import protein
from pybel_tools.selection import get_subgraph_by_data, get_subgraph_by_induction
from tests.constants import ExampleNetworkMixin, n

HGNC = 'HGNC'
GOBP = 'GOBP'
CHEBI = 'CHEBI'

g1 = GENE, HGNC, '1'
r1 = RNA, HGNC, '1'
p1 = PROTEIN, HGNC, '1'

g2 = GENE, HGNC, '2'
r2 = RNA, HGNC, '2'
p2 = PROTEIN, HGNC, '2'

g3 = GENE, HGNC, '3'
r3 = RNA, HGNC, '3'
p3 = PROTEIN, HGNC, '3'

g4 = GENE, HGNC, '4'
m4 = MIRNA, HGNC, '4'

a5 = ABUNDANCE, CHEBI, '5'
p5 = PATHOLOGY, GOBP, '5'


class TestInduction(unittest.TestCase):
    def test_induction_missing(self):
        """Checks that induction method returns none when no nodes present"""
        node1 = protein(namespace=n(), name=n())
        node2 = protein(namespace=n(), name=n())
        node3 = protein(namespace=n(), name=n())
        node4 = protein(namespace=n(), name=n())
        graph = BELGraph()
        graph.add_qualified_edge(node1, node2, relation=n(), citation=n(), evidence=n())
        graph.add_qualified_edge(node1, node4, relation=n(), citation=n(), evidence=n())

        res = get_subgraph_by_induction(graph, [node1.as_tuple()])
        self.assertIsNotNone(res)
        self.assertIsInstance(res, BELGraph)
        self.assertEqual(1, res.number_of_nodes())
        self.assertEqual(0, res.number_of_edges())

        res = get_subgraph_by_induction(graph, [node1.as_tuple(), node2.as_tuple()])
        self.assertIsNotNone(res)
        self.assertIsInstance(res, BELGraph)
        self.assertEqual(2, res.number_of_nodes())
        self.assertEqual(1, res.number_of_edges())

        res = get_subgraph_by_induction(graph, [node3.as_tuple()])
        self.assertIsNone(res, msg='Should return none since node3 is not in graph')

        self.assertEqual(3, graph.number_of_nodes(), msg='Original graph nodes should not change')
        self.assertEqual(2, graph.number_of_edges(), msg='Original graph edges should not change')


class TestGetSubgraphByData(ExampleNetworkMixin):
    def test_annotation_filter(self):
        test_network = self.network2  # Defined in test.constants.TestNetworks

        filtered_network = get_subgraph_by_data(test_network, {ANNOTATIONS: {'Annotation': 'foo2'}})

        self.assertEqual(2, filtered_network.number_of_nodes())
        self.assertEqual(1, filtered_network.number_of_edges())
        self.assertIn((GENE, 'HGNC', 'f'), filtered_network)
        self.assertIn((PROTEIN, 'HGNC', 'e'), filtered_network)

    def test_evidence_filter(self):
        test_network = self.graph_1  # Defined in test.constants.TestNetworks

        filtered_network = get_subgraph_by_data(test_network, {EVIDENCE: ['Evidence 1', 'Evidence 2']})

        self.assertEqual(3, filtered_network.number_of_nodes())
        self.assertEqual(2, filtered_network.number_of_edges())
        self.assertIn((PROTEIN, 'HGNC', 'a'), filtered_network)
        self.assertIn((PROTEIN, 'HGNC', 'b'), filtered_network)
        self.assertIn((RNA, 'HGNC', 'd'), filtered_network)

    def test_citation_filter(self):
        test_network = self.graph_1  # Defined in test.constants.TestNetworks

        filtered_network = get_subgraph_by_data(
            test_network,
            {
                CITATION: {
                    CITATION_TYPE: CITATION_TYPE_PUBMED,
                    CITATION_REFERENCE: ['2', '1'],
                },
            }
        )

        self.assertEqual(3, filtered_network.number_of_nodes())
        self.assertEqual(2, filtered_network.number_of_edges())

        self.assertIn((PROTEIN, 'HGNC', 'a'), filtered_network)
        self.assertIn((PROTEIN, 'HGNC', 'b'), filtered_network)
        self.assertIn((RNA, 'HGNC', 'd'), filtered_network)
