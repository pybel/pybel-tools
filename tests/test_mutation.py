# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.constants import unqualified_edge_code
from pybel.dsl import protein
from pybel_tools.mutation import collapse_by_central_dogma, collapse_nodes
from pybel_tools.mutation.collapse import collapse_to_protein_interactions
from pybel_tools.mutation.inference import infer_central_dogmatic_transcriptions, infer_central_dogmatic_translations
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


class TestCollapseDownstream(unittest.TestCase):
    def test_collapse_1(self):
        graph = BELGraph()

        graph.add_simple_node(*p1)
        graph.add_simple_node(*p2)
        graph.add_simple_node(*p3)

        graph.add_edge(p1, p3, **{RELATION: INCREASES})
        graph.add_edge(p2, p3, **{RELATION: DIRECTLY_INCREASES})

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        d = {
            p1: {p2}
        }

        collapse_nodes(graph, d)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges(), msg=graph.edges(data=True, keys=True))

    def test_collapse_dogma_1(self):
        graph = BELGraph()

        graph.add_simple_node(*p1)
        graph.add_simple_node(*r1)

        graph.add_edge(r1, p1, key=unqualified_edge_code[TRANSLATED_TO], **{RELATION: TRANSLATED_TO})

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        collapse_by_central_dogma(graph)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

    def test_collapse_dogma_2(self):
        graph = BELGraph()

        graph.add_simple_node(*p1)
        graph.add_simple_node(*r1)
        graph.add_simple_node(*g1)

        graph.add_edge(r1, p1, key=unqualified_edge_code[TRANSLATED_TO], **{RELATION: TRANSLATED_TO})
        graph.add_edge(g1, r1, key=unqualified_edge_code[TRANSCRIBED_TO], **{RELATION: TRANSCRIBED_TO})

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        collapse_by_central_dogma(graph)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

    def test_collapse_dogma_3(self):
        graph = BELGraph()

        graph.add_simple_node(*r1)
        graph.add_simple_node(*g1)

        graph.add_edge(g1, r1, key=unqualified_edge_code[TRANSCRIBED_TO], **{RELATION: TRANSCRIBED_TO})

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        collapse_by_central_dogma(graph)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())


class TestInference(unittest.TestCase):
    def test_infer_1(self):
        graph = BELGraph()

        graph.add_simple_node(*p1)
        graph.add_simple_node(*g1)
        graph.add_simple_node(*p2)
        graph.add_simple_node(*g3)
        graph.add_simple_node(*m4)

        graph.add_edge(p1, p2, **{RELATION: INCREASES})
        graph.add_edge(g1, g3, **{RELATION: POSITIVE_CORRELATION})

        self.assertEqual(5, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        infer_central_dogmatic_translations(graph)

        self.assertEqual(7, graph.number_of_nodes())
        self.assertEqual(4, graph.number_of_edges())
        self.assertIn(r1, graph)
        self.assertIn(r2, graph)

        infer_central_dogmatic_transcriptions(graph)

        self.assertEqual(9, graph.number_of_nodes())
        self.assertEqual(7, graph.number_of_edges())
        self.assertIn(g1, graph)
        self.assertIn(g2, graph)
        self.assertIn(g3, graph)
        self.assertIn(g4, graph)

        collapse_by_central_dogma(graph)

        self.assertEqual(4, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        self.assertTrue(graph.has_edge(p1, g3))
        self.assertTrue(graph.has_edge(p1, p2))


class TestGetSubgraphByData(ExampleNetworkMixin):
    def test_annotation_filter(self):
        test_network = self.network2  # Defined in test.constants.TestNetworks

        filtered_network = get_subgraph_by_data(test_network, {ANNOTATIONS: {'Annotation': 'foo2'}})

        self.assertEqual(2, filtered_network.number_of_nodes())
        self.assertEqual(1, filtered_network.number_of_edges())
        self.assertIn(('Gene', 'HGNC', 'f'), filtered_network)
        self.assertIn(('Protein', 'HGNC', 'e'), filtered_network)

    def test_evidence_filter(self):
        test_network = self.network1  # Defined in test.constants.TestNetworks

        filtered_network = get_subgraph_by_data(test_network, {EVIDENCE: ['Evidence 1', 'Evidence 2']})

        self.assertEqual(3, filtered_network.number_of_nodes())
        self.assertEqual(2, filtered_network.number_of_edges())
        self.assertIn(('Protein', 'HGNC', 'a'), filtered_network)
        self.assertIn(('Protein', 'HGNC', 'b'), filtered_network)
        self.assertIn(('RNA', 'HGNC', 'd'), filtered_network)

    def test_citation_filter(self):
        test_network = self.network1  # Defined in test.constants.TestNetworks

        filtered_network = get_subgraph_by_data(test_network,
                                                {
                                                    CITATION: {
                                                        CITATION_TYPE: CITATION_TYPE_PUBMED,
                                                        CITATION_REFERENCE: ['2', '1'],
                                                    },
                                                }
                                                )

        self.assertEqual(3, filtered_network.number_of_nodes())
        self.assertEqual(2, filtered_network.number_of_edges())

        self.assertIn(('Protein', 'HGNC', 'a'), filtered_network)
        self.assertIn(('Protein', 'HGNC', 'b'), filtered_network)
        self.assertIn(('RNA', 'HGNC', 'd'), filtered_network)


class TestCollapseProteinInteractions(unittest.TestCase):
    def test_protein_interaction_1(self):
        graph = BELGraph()

        graph.add_simple_node(*p1)
        graph.add_simple_node(*p2)
        graph.add_simple_node(*a5)
        graph.add_simple_node(*p5)

        graph.add_edge(p1, p2, **{RELATION: POSITIVE_CORRELATION})
        graph.add_edge(p1, p2, **{RELATION: INCREASES})
        graph.add_edge(a5, p5, **{RELATION: DIRECTLY_INCREASES})
        graph.add_edge(p1, a5, **{RELATION: DECREASES})

        collapsed_graph = collapse_to_protein_interactions(graph)

        self.assertEqual(2, collapsed_graph.number_of_nodes())
        self.assertEqual(2, collapsed_graph.number_of_edges())
        self.assertIn(('Gene', 'HGNC', '1'), collapsed_graph)
        self.assertIn(('Gene', 'HGNC', '2'), collapsed_graph)

    def test_protein_interaction_2(self):
        graph = BELGraph()

        graph.add_simple_node(*p1)
        graph.add_simple_node(*p2)
        graph.add_simple_node(*a5)
        graph.add_simple_node(*p5)

        graph.add_edge(p1, p2, **{RELATION: POSITIVE_CORRELATION})
        graph.add_edge(p1, p2, **{RELATION: ASSOCIATION})
        graph.add_edge(a5, p5, **{RELATION: DIRECTLY_INCREASES})
        graph.add_edge(p1, a5, **{RELATION: DECREASES})

        collapsed_graph = collapse_to_protein_interactions(graph)

        self.assertEqual(2, collapsed_graph.number_of_nodes())
        self.assertEqual(1, collapsed_graph.number_of_edges())
        self.assertIn(('Gene', 'HGNC', '1'), collapsed_graph)
        self.assertIn(('Gene', 'HGNC', '2'), collapsed_graph)
