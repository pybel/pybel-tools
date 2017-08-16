import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel_tools.mutation.inference import infer_missing_two_way_edges


class TestMutationInference(unittest.TestCase):
    def test_infer_missing_two_way_edges(self):
        g = BELGraph()

        a = PROTEIN, 'HGNC', 'A'
        b = PROTEIN, 'HGNC', 'B'
        c = PROTEIN, 'HGNC', 'C'
        d = PROTEIN, 'HGNC', 'D'

        g.add_simple_node(*a)
        g.add_simple_node(*b)
        g.add_simple_node(*c)
        g.add_simple_node(*d)

        g.add_edge(a, b, **{RELATION: POSITIVE_CORRELATION})
        g.add_edge(a, c, **{RELATION: POSITIVE_CORRELATION})
        g.add_edge(c, b, **{RELATION: NEGATIVE_CORRELATION})
        g.add_edge(a, d, **{RELATION: INCREASES})

        infer_missing_two_way_edges(g)

        self.assertTrue(g.has_edge(b, a))
        self.assertTrue(g.has_edge(c, a))
        self.assertTrue(g.has_edge(c, b))
        self.assertFalse(g.has_edge(d, a))
