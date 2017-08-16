import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel_tools.analysis.stability import *
from pybel_tools.mutation.inference import infer_missing_two_way_edges


class TestUnstableTriplets(unittest.TestCase):
    def test_separate_unstable(self):
        graph = BELGraph()

        a = PROTEIN, 'HGNC', 'A'
        b = PROTEIN, 'HGNC', 'B'
        c = PROTEIN, 'HGNC', 'C'
        d = PROTEIN, 'HGNC', 'D'

        graph.add_simple_node(*a)
        graph.add_simple_node(*b)
        graph.add_simple_node(*c)
        graph.add_simple_node(*d)

        graph.add_edge(a, b, **{RELATION: POSITIVE_CORRELATION})
        graph.add_edge(a, c, **{RELATION: POSITIVE_CORRELATION})
        graph.add_edge(c, b, **{RELATION: NEGATIVE_CORRELATION})

        infer_missing_two_way_edges(graph)

        cg = get_correlation_graph(graph)

        self.assertIn(a, cg)
        self.assertIn(b, cg)
        self.assertIn(c, cg)
        self.assertTrue(cg.has_edge(a, b))
        self.assertTrue(cg.has_edge(a, c))
        self.assertTrue(cg.has_edge(b, c))
        self.assertIn(POSITIVE_CORRELATION, cg.edge[a][b])
        self.assertIn(POSITIVE_CORRELATION, cg.edge[a][c])
        self.assertIn(NEGATIVE_CORRELATION, cg.edge[c][b])

        triangles = tuple(get_correlation_triangles(cg))

        self.assertEqual(1, len(triangles))
        self.assertEqual((a, b, c), triangles[0])

        result = tuple(get_separate_unstable_correlation_triples(graph))
        self.assertEqual(1, len(result))
        self.assertEqual((a, b, c), result[0])

    def test_mutually_unstable(self):
        g = BELGraph()

        a = PROTEIN, 'HGNC', 'A'
        b = PROTEIN, 'HGNC', 'B'
        c = PROTEIN, 'HGNC', 'C'
        d = PROTEIN, 'HGNC', 'D'

        g.add_simple_node(*a)
        g.add_simple_node(*b)
        g.add_simple_node(*c)
        g.add_simple_node(*d)

        g.add_edge(a, b, **{RELATION: NEGATIVE_CORRELATION})
        g.add_edge(a, c, **{RELATION: NEGATIVE_CORRELATION})
        g.add_edge(c, b, **{RELATION: NEGATIVE_CORRELATION})
        g.add_edge(c, d, **{RELATION: POSITIVE_CORRELATION})

        infer_missing_two_way_edges(g)

        result = tuple(get_mutually_unstable_correlation_triples(g))

        self.assertEqual(1, len(result))
        self.assertEqual((a, b, c), result[0])

    def test_jens_alpha(self):
        g = BELGraph()

        a = PROTEIN, 'HGNC', 'A'
        b = PROTEIN, 'HGNC', 'B'
        c = PROTEIN, 'HGNC', 'C'
        d = PROTEIN, 'HGNC', 'D'
        e = PROTEIN, 'HGNC', 'e'

        g.add_simple_node(*a)
        g.add_simple_node(*b)
        g.add_simple_node(*c)
        g.add_simple_node(*d)
        g.add_simple_node(*e)

        g.add_edge(a, b, **{RELATION: INCREASES})
        g.add_edge(a, c, **{RELATION: DECREASES})
        g.add_edge(c, b, **{RELATION: NEGATIVE_CORRELATION})
        g.add_edge(e, c, **{RELATION: POSITIVE_CORRELATION})
        g.add_edge(e, b, **{RELATION: POSITIVE_CORRELATION})
