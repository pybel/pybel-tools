# -*- coding: utf-8 -*-

from pybel.constants import (
    ABUNDANCE, ANNOTATIONS, CITATION, CITATION_REFERENCE, CITATION_TYPE, CITATION_TYPE_PUBMED, EVIDENCE, GENE, MIRNA,
    PATHOLOGY, PROTEIN, RNA,
)
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


class TestGetSubgraphByData(ExampleNetworkMixin):
    """Test building sub-graphs by data filters."""

    def test_annotation_filter(self):
        """Test building a sub-graph by an annotation filter."""
        test_network = self.network2  # Defined in test.constants.TestNetworks

        filtered_network = get_subgraph_by_data(test_network, {ANNOTATIONS: {'Annotation': 'foo2'}})

        self.assertEqual(2, filtered_network.number_of_nodes())
        self.assertEqual(1, filtered_network.number_of_edges())
        self.assertIn((GENE, 'HGNC', 'f'), filtered_network)
        self.assertIn((PROTEIN, 'HGNC', 'e'), filtered_network)

    def test_evidence_filter(self):
        """Test building a sub-graph by an evidence filter."""
        test_network = self.graph_1  # Defined in test.constants.TestNetworks

        filtered_network = get_subgraph_by_data(test_network, {EVIDENCE: ['Evidence 1', 'Evidence 2']})

        self.assertEqual(3, filtered_network.number_of_nodes())
        self.assertEqual(2, filtered_network.number_of_edges())
        self.assertIn((PROTEIN, 'HGNC', 'a'), filtered_network)
        self.assertIn((PROTEIN, 'HGNC', 'b'), filtered_network)
        self.assertIn((RNA, 'HGNC', 'd'), filtered_network)

    def test_citation_filter(self):
        """Test building a sub-graph by a citation filter."""
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
