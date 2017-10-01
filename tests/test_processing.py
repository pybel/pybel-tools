# -*- coding: utf-8 -*-

import os

import pybel
from pybel.constants import *
from pybel_tools.mutation import prune_central_dogma, infer_missing_inverse_edge, infer_central_dogma
from tests.constants import ManagerMixin


class TestProcessing(ManagerMixin):
    def setUp(self):
        super(TestProcessing, self).setUp()

        if 'PYBEL_BASE' in os.environ:
            test_bel_simple_path = os.path.join(os.environ['PYBEL_BASE'], 'tests', 'bel', 'test_bel.bel')
            self.graph = pybel.from_path(test_bel_simple_path, manager=self.manager)
        else:
            test_bel_simple_url = 'https://raw.githubusercontent.com/pybel/pybel/develop/tests/bel/test_bel.bel'
            self.graph = pybel.from_url(test_bel_simple_url, manager=self.manager)

        infer_central_dogma(self.graph)

        n1 = GENE, 'HGNC', 'AKT1'
        n2 = RNA, 'HGNC', 'EGFR'

        n3 = GENE, 'HGNC', 'DUMMY1'
        self.graph.add_simple_node(*n3)

        n4 = GENE, 'HGNC', 'DUMMY2'
        self.graph.add_simple_node(*n4)

        self.graph.add_edge(n1, n3, **{RELATION: INCREASES})
        self.graph.add_edge(n2, n4, **{RELATION: INCREASES})

    def test_prune(self):
        self.assertEqual(14, self.graph.number_of_nodes())
        self.assertEqual(16, self.graph.number_of_edges())

        prune_central_dogma(self.graph)
        self.assertEqual(9, self.graph.number_of_nodes())

    def test_infer(self):
        self.assertEqual(14, self.graph.number_of_nodes())
        self.assertEqual(16, self.graph.number_of_edges())

        infer_missing_inverse_edge(self.graph, TRANSLATED_TO)
        self.assertEqual(20, self.graph.number_of_edges())
