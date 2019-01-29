# -*- coding: utf-8 -*-

import unittest

from pybel.examples.sialic_acid_example import cd33, cd33_phosphorylated, sialic_acid_graph
from pybel.examples.various_example import (
    complex_example, composite_example, glycolisis_step_1, hk1, single_complex_graph, single_composite_graph,
    single_reaction_graph,
)
from pybel_tools.mutation import enrich_complexes, enrich_composites, enrich_reactions, enrich_variants


class TestEnrich(unittest.TestCase):
    def setUp(self):
        self.sialic_acid_graph = sialic_acid_graph.copy()
        self.single_reaction_graph = single_reaction_graph.copy()
        self.single_composite_graph = single_composite_graph.copy()
        self.single_complex_graph = single_complex_graph.copy()

    def test_enrich_reactions(self):
        """Test enrich enrichment."""
        self.assertEqual(7, self.single_reaction_graph.number_of_nodes())
        self.assertEqual(7, self.single_reaction_graph.number_of_edges())

        self.single_reaction_graph.remove_edge(glycolisis_step_1, hk1)

        self.assertEqual(7, self.single_reaction_graph.number_of_nodes())
        self.assertEqual(6, self.single_reaction_graph.number_of_edges())

        enrich_reactions(self.single_reaction_graph)

        self.assertEqual(7, self.single_reaction_graph.number_of_nodes(), msg='Enrich reactions did not work')
        self.assertEqual(7, self.single_reaction_graph.number_of_edges(), msg='Enrich reactions did not work')

    def test_enrich_composites(self):
        """Test composite enrichment."""
        self.assertEqual(4, self.single_composite_graph.number_of_nodes())
        self.assertEqual(3, self.single_composite_graph.number_of_edges())

        self.single_composite_graph.remove_edge(composite_example, hk1)

        self.assertEqual(4, self.single_composite_graph.number_of_nodes())
        self.assertEqual(2, self.single_composite_graph.number_of_edges())

        enrich_composites(self.single_composite_graph)

        self.assertEqual(4, self.single_composite_graph.number_of_nodes(), msg='Enrich composites did not work')
        self.assertEqual(3, self.single_composite_graph.number_of_edges(), msg='Enrich composites did not work')

    def test_enrich_complexes(self):
        """Test composite enrichment."""
        self.assertEqual(4, self.single_complex_graph.number_of_nodes())
        self.assertEqual(3, self.single_complex_graph.number_of_edges())

        self.single_complex_graph.remove_edge(complex_example, hk1)

        self.assertEqual(4, self.single_complex_graph.number_of_nodes())
        self.assertEqual(2, self.single_complex_graph.number_of_edges())

        enrich_complexes(self.single_complex_graph)

        self.assertEqual(4, self.single_complex_graph.number_of_nodes(), msg='Enrich complexes did not work')
        self.assertEqual(3, self.single_complex_graph.number_of_edges(), msg='Enrich complexes did not work')

    def test_enrich_variants(self):
        """Test enrich variants."""
        self.sialic_acid_graph.remove_node(cd33)

        self.assertNotIn(cd33, self.sialic_acid_graph)
        self.assertIn(cd33_phosphorylated, self.sialic_acid_graph)

        enrich_variants(self.sialic_acid_graph)

        self.assertIn(cd33_phosphorylated, self.sialic_acid_graph)
        self.assertIn(cd33, self.sialic_acid_graph, msg='Enrich variants did not work')
