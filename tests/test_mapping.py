import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel_tools.selection.group_nodes import get_mapped


class TestMapping(unittest.TestCase):
    def test_variants_mapping(self):
        g = BELGraph()

        g.add_node('APP', attr_dict={NAMESPACE: 'HGNC', NAME: 'APP'})
        g.add_node('APP Fragment')
        g.add_edge('APP', 'APP Fragment', **{RELATION: HAS_VARIANT})

        mapped_nodes = get_mapped(g, 'HGNC', {'APP'})

        self.assertEqual(1, len(mapped_nodes))
        self.assertIn('APP', mapped_nodes)
        self.assertEqual({'APP Fragment'}, mapped_nodes['APP'])

    def test_complexes_composites_mapping(self):
        g = BELGraph()

        g.add_node('complex(p(HGNC:CCL2), p(HGNC:CCR2))')
        g.add_node('CCL2', attr_dict={NAMESPACE: 'HGNC', NAME: 'CCL2'})
        g.add_node('CCR2', attr_dict={NAMESPACE: 'HGNC', NAME: 'CCR2'})

        g.add_node('chemokine protein family')

        g.add_edge('chemokine protein family', 'CCL2', **{RELATION: HAS_MEMBER})
        g.add_edge('chemokine protein family', 'CCR2', **{RELATION: HAS_MEMBER})

        g.add_edge('complex(p(HGNC:CCL2), p(HGNC:CCR2))', 'CCL2', **{RELATION: HAS_COMPONENT})
        g.add_edge('complex(p(HGNC:CCL2), p(HGNC:CCR2))', 'CCR2', **{RELATION: HAS_COMPONENT})

        mapped_nodes = get_mapped(g, 'HGNC', {'CCL2', 'CCR2'})

        self.assertEqual(2, len(mapped_nodes))
        self.assertIn('CCL2', mapped_nodes)
        self.assertIn('CCR2', mapped_nodes)

        self.assertEqual({'complex(p(HGNC:CCL2), p(HGNC:CCR2))', 'chemokine protein family'}, mapped_nodes['CCR2'])
        self.assertEqual({'complex(p(HGNC:CCL2), p(HGNC:CCR2))', 'chemokine protein family'}, mapped_nodes['CCR2'])


    def test_orthologus_mapping(self):
        g = BELGraph()

        g.add_node('Ccl2', attr_dict={NAMESPACE: 'MGI', NAME: 'Ccl2'})
        g.add_node('CCL2', attr_dict={NAMESPACE: 'HGNC', NAME: 'CCL2'})

        g.add_edge('CCL2', 'Ccl2', **{RELATION: ORTHOLOGOUS})

        mapped_nodes = get_mapped(g, 'HGNC', {'CCL2'})

        self.assertEqual(1, len(mapped_nodes))
        self.assertIn('CCL2', mapped_nodes)

        self.assertEqual({'Ccl2'}, mapped_nodes['CCL2'])