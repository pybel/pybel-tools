# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.dsl import ComplexAbundance, Fragment, Protein
from pybel.dsl.namespaces import hgnc
from pybel_tools.selection.group_nodes import get_mapped_nodes

ccl2 = hgnc(name='CCL2')
ccr2 = hgnc(name='CCR2')
ccl2_mgi = Protein('MGI', 'Ccl2')
ccl2_ccr2_complex = ComplexAbundance([ccl2, ccr2])
chemokine_family = Protein('FPLX', 'chemokine protein family')

HGNC = 'hgnc'


class TestMapping(unittest.TestCase):
    def test_variants_mapping(self):
        graph = BELGraph()

        app = Protein(HGNC, 'APP')
        app_fragment = app.with_variants(Fragment('1_49'))
        graph.add_node_from_data(app_fragment)

        mapped_nodes = get_mapped_nodes(graph, HGNC, {'APP'})

        self.assertEqual(1, len(mapped_nodes))
        self.assertIn(app, mapped_nodes)
        self.assertEqual({app_fragment}, mapped_nodes[app])

    def test_complexes_composites_mapping(self):
        g = BELGraph()
        g.add_is_a(ccl2, chemokine_family)
        g.add_is_a(ccr2, chemokine_family)
        g.add_part_of(ccl2, ccl2_ccr2_complex)
        g.add_part_of(ccr2, ccl2_ccr2_complex)

        mapped_nodes = get_mapped_nodes(g, 'HGNC', {ccl2.name, ccr2.name})

        self.assertEqual(2, len(mapped_nodes))
        self.assertIn(ccl2, mapped_nodes)
        self.assertIn(ccr2, mapped_nodes)

        self.assertEqual({ccl2_ccr2_complex, chemokine_family}, mapped_nodes[ccl2])
        self.assertEqual({ccl2_ccr2_complex, chemokine_family}, mapped_nodes[ccr2])

    def test_orthologus_mapping(self):
        g = BELGraph()

        g.add_node_from_data(ccl2)
        g.add_node_from_data(ccl2_mgi)
        g.add_orthology(ccl2, ccl2_mgi)

        mapped_nodes = get_mapped_nodes(g, 'HGNC', {'CCL2'})

        self.assertEqual(1, len(mapped_nodes))
        self.assertIn(ccl2, mapped_nodes)

        self.assertEqual({ccl2_mgi}, mapped_nodes[ccl2])
