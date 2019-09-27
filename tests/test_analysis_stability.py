import unittest

from pybel import BELGraph
from pybel.constants import DECREASES, INCREASES, NEGATIVE_CORRELATION, POSITIVE_CORRELATION, RELATION
from pybel.dsl import Protein
from pybel_tools.mutation.inference import infer_missing_two_way_edges
from pybel_tools.summary import (
    get_correlation_graph, get_correlation_triangles, get_mutually_unstable_correlation_triples,
    get_separate_unstable_correlation_triples,
)


class TestUnstableTriplets(unittest.TestCase):
    def test_separate_unstable(self):
        graph = BELGraph()

        a_data = Protein('HGNC', 'A')
        b_data = Protein('HGNC', 'B')
        c_data = Protein('HGNC', 'C')
        d_data = Protein('HGNC', 'D')

        a = graph.add_node_from_data(a_data)
        b = graph.add_node_from_data(b_data)
        c = graph.add_node_from_data(c_data)
        d = graph.add_node_from_data(d_data)

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
        self.assertIn(POSITIVE_CORRELATION, cg[a][b])
        self.assertIn(POSITIVE_CORRELATION, cg[a][c])
        self.assertIn(NEGATIVE_CORRELATION, cg[c][b])

        triangles = tuple(get_correlation_triangles(cg))

        self.assertEqual(1, len(triangles))
        self.assertEqual((a, b, c), triangles[0])

        result = tuple(get_separate_unstable_correlation_triples(graph))
        self.assertEqual(1, len(result))
        self.assertEqual((a, b, c), result[0])

    def test_mutually_unstable(self):
        graph = BELGraph()

        a_data = Protein('HGNC', 'A')
        b_data = Protein('HGNC', 'B')
        c_data = Protein('HGNC', 'C')
        d_data = Protein('HGNC', 'D')

        a = graph.add_node_from_data(a_data)
        b = graph.add_node_from_data(b_data)
        c = graph.add_node_from_data(c_data)
        d = graph.add_node_from_data(d_data)

        graph.add_edge(a, b, **{RELATION: NEGATIVE_CORRELATION})
        graph.add_edge(a, c, **{RELATION: NEGATIVE_CORRELATION})
        graph.add_edge(c, b, **{RELATION: NEGATIVE_CORRELATION})
        graph.add_edge(c, d, **{RELATION: POSITIVE_CORRELATION})

        infer_missing_two_way_edges(graph)

        result = tuple(get_mutually_unstable_correlation_triples(graph))

        self.assertEqual(1, len(result))
        self.assertEqual((a, b, c), result[0])

    def test_jens_alpha(self):
        graph = BELGraph()

        a_data = Protein('HGNC', 'A')
        b_data = Protein('HGNC', 'B')
        c_data = Protein('HGNC', 'C')
        d_data = Protein('HGNC', 'D')
        e_data = Protein('HGNC', 'e')

        a = graph.add_node_from_data(a_data)
        b = graph.add_node_from_data(b_data)
        c = graph.add_node_from_data(c_data)
        d = graph.add_node_from_data(d_data)
        e = graph.add_node_from_data(e_data)

        graph.add_edge(a, b, **{RELATION: INCREASES})
        graph.add_edge(a, c, **{RELATION: DECREASES})
        graph.add_edge(c, b, **{RELATION: NEGATIVE_CORRELATION})
        graph.add_edge(e, c, **{RELATION: POSITIVE_CORRELATION})
        graph.add_edge(e, b, **{RELATION: POSITIVE_CORRELATION})
