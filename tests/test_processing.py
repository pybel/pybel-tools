# -*- coding: utf-8 -*-

import unittest

from pybel.dsl import gene, rna
from pybel.examples import sialic_acid_graph
from pybel_tools.mutation import infer_central_dogma, prune_central_dogma

trem2_gene = gene(namespace='HGNC', name='TREM2')
trem2_rna = rna(namespace='HGNC', name='TREM2')


class TestProcessing(unittest.TestCase):
    def test_infer_on_sialic_acid_example(self):
        sgc = sialic_acid_graph.copy()

        self.assertFalse(sgc.has_node_with_data(trem2_gene))
        self.assertFalse(sgc.has_node_with_data(trem2_rna))

        infer_central_dogma(sgc)

        self.assertTrue(sgc.has_node_with_data(trem2_gene))
        self.assertTrue(sgc.has_node_with_data(trem2_rna))

        prune_central_dogma(sgc)

        self.assertFalse(sgc.has_node_with_data(trem2_gene))
        self.assertFalse(sgc.has_node_with_data(trem2_rna))
