# -*- coding: utf-8 -*-

import os
import tempfile
import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.manager import CacheManager

HGNC = 'HGNC'

dir_path = os.path.dirname(os.path.realpath(__file__))
resources_path = os.path.join(dir_path, 'resources')

rgd_orthologs_path = os.path.join(resources_path, 'RGD_ORTHOLOGS.txt')


class ManagerMixin(unittest.TestCase):
    def setUp(self):
        super(ManagerMixin, self).setUp()

        self.db_fd, self.db_file = tempfile.mkstemp()

        self.connection = 'sqlite:///' + self.db_file
        self.manager = CacheManager(connection=self.connection)

    def tearDown(self):
        os.close(self.db_fd)
        os.unlink(self.db_file)


protein_a = PROTEIN, HGNC, 'a'
protein_b = PROTEIN, HGNC, 'b'
gene_c = GENE, HGNC, 'c'
rna_d = RNA, HGNC, 'd'
protein_e = PROTEIN, HGNC, 'e'
gene_f = GENE, HGNC, 'f'


def make_graph_1():
    graph1 = BELGraph(**{
        GRAPH_METADATA: {
            METADATA_VERSION: '1.1.0',
            METADATA_NAME: 'network_test',
            METADATA_DESCRIPTION: 'network test',
            METADATA_AUTHORS: 'Fraunhofer SCAI',
            METADATA_CONTACT: 'test@scai.fraunhofer.de',
        }
    })

    graph1.add_simple_node(*protein_a)
    graph1.add_simple_node(*protein_b)
    graph1.add_simple_node(*gene_c)
    graph1.add_simple_node(*rna_d)

    graph1.add_edge(protein_a, protein_b, attr_dict={
        RELATION: INCREASES,
        CITATION: {
            CITATION_TYPE: CITATION_TYPE_PUBMED,
            CITATION_REFERENCE: '1',
        },
        EVIDENCE: 'Evidence 1',
        ANNOTATIONS: {'Annotation': 'foo'}
    })

    graph1.add_edge(rna_d, protein_a, attr_dict={
        RELATION: INCREASES,
        CITATION: {
            CITATION_TYPE: CITATION_TYPE_PUBMED,
            CITATION_REFERENCE: '2',
        },
        EVIDENCE: 'Evidence 2',
        ANNOTATIONS: {'Annotation': 'foo'}
    })

    graph1.add_edge(gene_c, protein_b, attr_dict={
        RELATION: DECREASES,
        CITATION: {
            CITATION_TYPE: CITATION_TYPE_PUBMED,
            CITATION_REFERENCE: '3',
        },
        EVIDENCE: 'Evidence 3',
        ANNOTATIONS: {'Annotation': 'foo'}
    })

    return graph1


def make_graph_2():
    graph2 = BELGraph(**{
        GRAPH_METADATA: {
            METADATA_VERSION: '1.0.0',
            METADATA_NAME: 'network_test',
            METADATA_DESCRIPTION: 'network test',
            METADATA_AUTHORS: 'Fraunhofer SCAI',
            METADATA_CONTACT: 'test@scai.fraunhofer.de',
        }
    })

    graph2.add_simple_node(*gene_f)
    graph2.add_simple_node(*protein_e)
    graph2.add_simple_node(*protein_b)

    graph2.add_edge(protein_e, protein_b, attr_dict={
        RELATION: INCREASES,
        CITATION: {
            CITATION_TYPE: CITATION_TYPE_PUBMED,
            CITATION_REFERENCE: '1',
        },
        EVIDENCE: 'Evidence 1',
        ANNOTATIONS: {'Annotation': 'foo'}
    })

    graph2.add_edge(gene_f, protein_e, attr_dict={
        RELATION: INCREASES,
        CITATION: {
            CITATION_TYPE: CITATION_TYPE_PUBMED,
            CITATION_REFERENCE: '2',
        },
        EVIDENCE: 'Evidence 2',
        ANNOTATIONS: {'Annotation': 'foo2'}
    })

    return graph2


class ExampleNetworkMixin(unittest.TestCase):
    def setUp(self):
        super(ExampleNetworkMixin, self).setUp()
        self.network1 = make_graph_1()
        self.network2 = make_graph_2()
