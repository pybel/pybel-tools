# -*- coding: utf-8 -*-


import logging
from pybel.constants import *
from pybel_tools.analysis.sst import run_cna, rank_causalr_hypothesis
from tests.constants import ExampleNetworkMixin, ManagerMixin
from pybel.dsl import protein, rna, gene

HGNC = 'HGNC'

log = logging.getLogger(__name__)
log.setLevel(10)


class SstTest(ExampleNetworkMixin):
    def test_cna(self):
        test_network = self.network3  # Defined in test.constants.TestNetworks

        a_as_root = run_cna(graph=test_network, root=protein( HGNC, 'a'), targets=[gene( HGNC, 'c')])

        self.assertEqual(-1, a_as_root[0][2].value)  # A -| C

        d_as_root = run_cna(graph=test_network, root=rna( HGNC, 'd'), targets=[gene( HGNC, 'c'), gene( HGNC, 'f')])

        self.assertEqual(-1, d_as_root[0][2].value)  # D -| C
        self.assertEqual(-1, d_as_root[1][2].value)  # A -| F

        e_as_root = run_cna(graph=test_network, root=protein( HGNC, 'e'),
                            targets=[gene( HGNC, 'c'), gene( HGNC, 'f')])

        self.assertEqual(1, e_as_root[0][2].value)  # E -> C
        self.assertEqual(1, e_as_root[1][2].value)  # E -> F

        failed_results = run_cna(graph=test_network, root=protein( HGNC, 'e'), targets=[protein( HGNC, 'g')])

        self.assertEqual(0, failed_results[0][2].value)  # E -> G

    def test_causalr_rank_hypothesis_1(self):
        test_network = self.network4  # Defined in test.constants.TestNetworks

        observed_regulation_test = {
            protein( HGNC, 'a'): 0,
            protein( HGNC, 'b'): 1,
            gene( HGNC, 'c'): -1,
            rna( HGNC, 'd'): -1,
            protein( HGNC, 'e'): -1,
            gene( HGNC, 'f'): 1,
            protein( HGNC, 'g'): 1,
            protein( HGNC, 'h'): 1,
            protein( HGNC, 'i'): 1,
            protein( HGNC, 'j'): 1
        }

        upregulated_hypothesis, downregulated_hypothesis = rank_causalr_hypothesis(
            graph=test_network,
            node_to_regulation=observed_regulation_test,
            regulator_node=protein( HGNC, 'a')
        )

        self.assertEqual(5, upregulated_hypothesis['score'])
        self.assertEqual(6, upregulated_hypothesis['correct'])
        self.assertEqual(1, upregulated_hypothesis['incorrect'])  # 1 gene( HGNC, 'f')
        self.assertEqual(1, upregulated_hypothesis['ambiguous'])  # 1 gene( HGNC, 'h')

    def test_causalr_rank_hypothesis_2(self):
        test_network = self.network4  # Defined in test.constants.TestNetworks

        observed_regulation_test = {
            protein( HGNC, 'a'): 0,
            protein( HGNC, 'b'): 1,
            gene( HGNC, 'c'): 1,
            rna( HGNC, 'd'): 1,
            protein( HGNC, 'e'): 1,
            gene( HGNC, 'f'): -1,
            protein( HGNC, 'g'): -1,
            protein( HGNC, 'h'): -1,
            protein( HGNC, 'i'): -1,
            protein( HGNC, 'j'): 1
        }

        upregulated_hypothesis, downregulated_hypothesis = rank_causalr_hypothesis(
            graph=test_network,
            node_to_regulation=observed_regulation_test,
            regulator_node=protein( HGNC, 'b')
        )

        self.assertEqual(4, downregulated_hypothesis['score'])
        self.assertEqual(5, downregulated_hypothesis['correct'])
        self.assertEqual(1, downregulated_hypothesis['incorrect'])
        self.assertEqual(1, downregulated_hypothesis['ambiguous'])  # 1 gene( HGNC, 'h')

        # Checking for scoring symmetry
        self.assertEqual(abs(downregulated_hypothesis['score']), abs(upregulated_hypothesis['score']))
        self.assertEqual(abs(downregulated_hypothesis['incorrect']), abs(upregulated_hypothesis['correct']))
        self.assertEqual(abs(downregulated_hypothesis['correct']), abs(upregulated_hypothesis['incorrect']))
        self.assertEqual(abs(downregulated_hypothesis['ambiguous']), abs(upregulated_hypothesis['ambiguous']))
