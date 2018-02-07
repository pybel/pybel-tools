# -*- coding: utf-8 -*-

import unittest

import pybel_tools.pipeline
from pybel import BELGraph
from pybel.constants import *
from pybel.constants import unqualified_edge_code
from pybel.examples import egf_example
from pybel_tools.mutation import collapse_by_central_dogma, collapse_nodes
from pybel_tools.mutation.bound import build_delete_node_by_hash, build_expand_node_neighborhood_by_hash
from pybel_tools.mutation.collapse import collapse_to_protein_interactions
from pybel_tools.mutation.inference import infer_central_dogmatic_transcriptions, infer_central_dogmatic_translations
from pybel_tools.pipeline import Pipeline
from pybel_tools.selection import get_subgraph_by_data
from tests.constants import ExampleNetworkMixin

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
                                                {CITATION: {
                                                    CITATION_TYPE: CITATION_TYPE_PUBMED,
                                                    CITATION_REFERENCE: ['2', '1'],
                                                }, }
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


class MockManager(object):
    def __init__(self, graph):
        """Instantiates a mock manager for testing bound mutation methods

        :param pybel.BELGraph graph: A BEL graph to index
        """
        self.graph = graph

        self.hash_to_tuple = {
            graph.hash_node(data): node
            for node, data in graph.nodes(data=True)
        }

    def get_node_tuple_by_hash(self, node_hash):
        return self.hash_to_tuple[node_hash]


class TestBoundMutation(unittest.TestCase):
    """Random test for mutation functions"""

    def test_bound_mutation(self):
        """Tests when a node is deleted then re-expanded"""
        graph = egf_example.egf_graph

        original_number_nodes = graph.number_of_nodes()
        original_number_edges = graph.number_of_edges()

        manager = MockManager(graph)

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
