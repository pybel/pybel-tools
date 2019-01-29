# -*- coding: utf-8 -*-

"""Tests for inference functions."""

import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.dsl import protein
from pybel_tools.mutation.inference import infer_missing_two_way_edges


class TestMutationInference(unittest.TestCase):
    def test_infer_missing_two_way_edges(self):
        graph = BELGraph()

        a_data = protein('HGNC', 'A')
        b_data = protein('HGNC', 'B')
        c_data = protein('HGNC', 'C')
        d_data = protein('HGNC', 'D')

        a = graph.add_node_from_data(a_data)
        b = graph.add_node_from_data(b_data)
        c = graph.add_node_from_data(c_data)
        d = graph.add_node_from_data(d_data)

        graph.add_edge(a, b, **{RELATION: POSITIVE_CORRELATION})
        graph.add_edge(a, c, **{RELATION: POSITIVE_CORRELATION})
        graph.add_edge(c, b, **{RELATION: NEGATIVE_CORRELATION})
        graph.add_edge(a, d, **{RELATION: INCREASES})

        infer_missing_two_way_edges(graph)

        self.assertTrue(graph.has_edge(b, a))
        self.assertTrue(graph.has_edge(c, a))
        self.assertTrue(graph.has_edge(c, b))
        self.assertFalse(graph.has_edge(d, a))
