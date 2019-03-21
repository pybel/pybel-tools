# -*- coding: utf-8 -*-

"""Tests for node utilities."""

import unittest

from pybel_tools.mutation.collapse import collapse_nodes_with_same_names
from pybel_tools.node_utils import flatten_list_abundance

from pybel import BELGraph
from pybel.dsl import ComplexAbundance as g, CompositeAbundance as c, Protein
from pybel.testing.utils import n


class TestNodeUtils(unittest.TestCase):
    """Test node utilities."""

    def test_flatten_complex(self):
        """Test flattening a nested complex."""
        p1, p2, p3 = (Protein('N', str(i + 1)) for i in range(3))

        pairs = [
            # Mainly complexes
            (g([p1, p2, p3]), g([p1, p2, p3])),  # no nesting
            (g([p1, p2, p3]), g([g([p1, p2]), p3])),  # one nesting
            (g([p1, p2, p3]), g([g([p1]), p2, p3])),  # one nesting
            (g([p1, p2, p3]), g([g([p1]), g([p2]), p3])),  # one nesting
            # Mainly composites
            (c([p1, p2, p3]), c([p1, p2, p3])),  # no nesting
            (c([p1, p2, p3]), c([c([p1, p2]), p3])),  # one nesting
            (c([p1, p2, p3]), c([c([p1]), p2, p3])),  # one nesting
            (c([p1, p2, p3]), c([c([p1]), c([p2]), p3])),  # one nesting
            # TODO: mixtures of composites and complexes?
        ]

        for expected, source in pairs:
            self.assertEqual(expected, flatten_list_abundance(source))

    def test_merge_nodes_by_name(self):
        graph = BELGraph()
        a_hgnc, a_entrez = [Protein(namespace, 'a') for namespace in ('HGNC', 'ncbigene')]
        b = Protein('ncbigene', 'b')
        graph.add_increases(a_hgnc, b, n(), n())
        graph.add_increases(a_entrez, b, n(), n())
        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())
        collapse_nodes_with_same_names(graph)
        self.assertEqual(2, graph.number_of_nodes(), msg='Wrong number remaining nodes')
        self.assertEqual(2, graph.number_of_edges(), msg=f'Wrong number remaining edges: {graph.edges()}')
        self.assertIn(a_entrez, graph)
        self.assertIn(b, graph)
        self.assertNotIn(a_hgnc, graph)
