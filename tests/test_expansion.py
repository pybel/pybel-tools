# -*- coding: utf-8 -*-

import logging
import unittest

from pybel.examples.sialic_acid_example import sialic_acid_graph, cd33, cd33_phosphorylated
from pybel.testing.mock_manager import MockQueryManager
from pybel_tools.mutation import enrich_variants

log = logging.getLogger(__name__)
log.setLevel(10)


class TestBoundMutation(unittest.TestCase):
    """Random test for mutation functions."""

    def setUp(self):
        self.graph = sialic_acid_graph.copy()
        self.original_number_nodes = self.graph.number_of_nodes()
        self.original_number_edges = self.graph.number_of_edges()

        self.manager = MockQueryManager()
        self.network_id = self.manager.insert_graph(self.graph).id

    def test_enrich_composite(self):
        """Test enrich grouping."""
        self.graph.remove_node(cd33)

        self.assertNotIn(cd33, self.graph)
        self.assertIn(cd33_phosphorylated, self.graph)

        enrich_variants(self.graph)

        self.assertIn(cd33_phosphorylated, self.graph)
        self.assertIn(cd33, self.graph, msg='Enrich variants did not work')

    def test_enrich_variants(self):
        """Test enrich variants."""
        self.graph.remove_node(cd33)

        self.assertNotIn(cd33, self.graph)
        self.assertIn(cd33_phosphorylated, self.graph)

        enrich_variants(self.graph)

        self.assertIn(cd33_phosphorylated, self.graph)
        self.assertIn(cd33, self.graph, msg='Enrich variants did not work')
