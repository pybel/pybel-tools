# -*- coding: utf-8 -*-

"""Tests for node utilities."""

import unittest

from pybel.dsl import ComplexAbundance as g, CompositeAbundance as c, Protein
from pybel_tools.node_utils import flatten_list_abundance


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
