# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.dsl import complex_abundance, fragment, protein
from pybel.dsl.namespaces import hgnc
from pybel_tools.selection.group_nodes import get_mapped_nodes

ccl2 = hgnc('CCL2')
ccr2 = hgnc('CCR2')
ccl2_mgi = protein('MGI', 'Ccl2')
ccl2_ccr2_complex = complex_abundance([ccl2, ccr2])
chemokine_family = protein('FPLX', 'chemokine protein family')

HGNC = 'HGNC'


class TestMapping(unittest.TestCase):
    def test_variants_mapping(self):
        graph = BELGraph()

        app = protein(HGNC, 'APP')
        app_fragment = app.with_variants(fragment('1_49'))
        graph.add_node_from_data(app_fragment)

        mapped_nodes = get_mapped_nodes(graph, HGNC, {'APP'})

        self.assertEqual(1, len(mapped_nodes))
        self.assertIn(app, mapped_nodes)
        self.assertEqual({app_fragment}, mapped_nodes[app])

    def test_complexes_composites_mapping(self):
        g = BELGraph()
        g.add_node_from_data(ccl2)
        g.add_node_from_data(ccr2)
        g.add_node_from_data(ccl2_ccr2_complex)
        g.add_node_from_data(chemokine_family)

        g.add_has_member(chemokine_family, ccl2)
        g.add_has_member(chemokine_family, ccr2)

        g.add_has_component(ccl2_ccr2_complex, ccl2)
        g.add_has_component(ccl2_ccr2_complex, ccr2)

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
