# -*- coding: utf-8 -*-

import os
import tempfile
import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.dsl import gene, protein, rna
from pybel.manager import Manager
from uuid import uuid4

HGNC = 'HGNC'

dir_path = os.path.dirname(os.path.realpath(__file__))


class ManagerMixin(unittest.TestCase):
    def setUp(self):
        super(ManagerMixin, self).setUp()

        self.db_fd, self.db_file = tempfile.mkstemp()

        self.connection = 'sqlite:///' + self.db_file
        self.manager = Manager(connection=self.connection)

    def tearDown(self):
        os.close(self.db_fd)
        os.unlink(self.db_file)


protein_a = protein(namespace=HGNC, name='a')
protein_b = protein(namespace=HGNC, name='b')
gene_c = gene(namespace=HGNC, name='c')
rna_d = rna(namespace=HGNC, name='d')
protein_e = protein(namespace=HGNC, name='e')
gene_f = gene(namespace=HGNC, name='f')
protein_g = protein(namespace=HGNC, name='g')
protein_h = protein(namespace=HGNC, name='h')
protein_i = protein(namespace=HGNC, name='i')
protein_j = protein(namespace=HGNC, name='j')

protein_a_tuple = PROTEIN, HGNC, 'a'
protein_b_tuple = PROTEIN, HGNC, 'b'
gene_c_tuple = GENE, HGNC, 'c'
rna_d_tuple = RNA, HGNC, 'd'
protein_e_tuple = PROTEIN, HGNC, 'e'
gene_f_tuple = GENE, HGNC, 'f'
protein_g_tuple = PROTEIN, HGNC, 'g'
protein_h_tuple = PROTEIN, HGNC, 'h'
protein_i_tuple = PROTEIN, HGNC, 'i'
protein_j_tuple = PROTEIN, HGNC, 'j'


def make_graph_1():
    graph = BELGraph(
        name='PyBEL Tools Example Network 1',
        version='1.1.0',
        description='Example Network for PyBEL Tools Tests',
        authors='Daniel Domingo-Fernández and Charles Tapley Hoyt',
        contact='charles.hoyt@scai.fraunhofer.de',
    )

    graph.add_node_from_data(protein_a)
    graph.add_node_from_data(protein_b)
    graph.add_node_from_data(gene_c)
    graph.add_node_from_data(rna_d)

    graph.add_edge(protein_a_tuple, protein_b_tuple, attr_dict={
        RELATION: INCREASES,
        CITATION: {
            CITATION_TYPE: CITATION_TYPE_PUBMED,
            CITATION_REFERENCE: '1',
        },
        EVIDENCE: 'Evidence 1',
        ANNOTATIONS: {'Annotation': 'foo'}
    })

    graph.add_edge(rna_d_tuple, protein_a_tuple, attr_dict={
        RELATION: INCREASES,
        CITATION: {
            CITATION_TYPE: CITATION_TYPE_PUBMED,
            CITATION_REFERENCE: '2',
        },
        EVIDENCE: 'Evidence 2',
        ANNOTATIONS: {'Annotation': 'foo'}
    })

    graph.add_edge(gene_c_tuple, protein_b_tuple, attr_dict={
        RELATION: DECREASES,
        CITATION: {
            CITATION_TYPE: CITATION_TYPE_PUBMED,
            CITATION_REFERENCE: '3',
        },
        EVIDENCE: 'Evidence 3',
        ANNOTATIONS: {'Annotation': 'foo'}
    })

    return graph


def make_graph_2():
    graph = BELGraph(
        name='PyBEL Tools Example Network 2',
        version='1.0.0',
        description='Example Network for PyBEL Tools Tests',
        authors='Daniel Domingo-Fernández and Charles Tapley Hoyt',
        contact='charles.hoyt@scai.fraunhofer.de',
    )

    graph.add_node_from_data(gene_f)
    graph.add_node_from_data(protein_e)
    graph.add_node_from_data(protein_b)

    graph.add_edge(protein_e_tuple, protein_b_tuple, attr_dict={
        RELATION: INCREASES,
        CITATION: {
            CITATION_TYPE: CITATION_TYPE_PUBMED,
            CITATION_REFERENCE: '1',
        },
        EVIDENCE: 'Evidence 1',
        ANNOTATIONS: {'Annotation': 'foo'}
    })

    graph.add_edge(gene_f_tuple, protein_e_tuple, attr_dict={
        RELATION: INCREASES,
        CITATION: {
            CITATION_TYPE: CITATION_TYPE_PUBMED,
            CITATION_REFERENCE: '2',
        },
        EVIDENCE: 'Evidence 2',
        ANNOTATIONS: {'Annotation': 'foo2'}
    })

    return graph


def make_graph_3():
    """
    A -> B -| C
    D -| F -> C
    C -| F
    C -- G

    """
    graph = BELGraph(
        name='PyBEL Tools Example Network 3',
        version='1.0.0',
        description='Example Network for PyBEL Tools Tests',
        authors='Daniel Domingo-Fernández and Charles Tapley Hoyt',
        contact='charles.hoyt@scai.fraunhofer.de',
    )

    graph.add_edge(protein_a_tuple, protein_b_tuple, attr_dict={
        RELATION: INCREASES,
    })

    graph.add_edge(protein_b_tuple, gene_c_tuple, attr_dict={
        RELATION: DECREASES,
    })

    graph.add_edge(rna_d_tuple, gene_f_tuple, attr_dict={
        RELATION: DECREASES,
    })

    graph.add_edge(protein_e_tuple, gene_f_tuple, attr_dict={
        RELATION: INCREASES,
    })

    graph.add_edge(gene_f_tuple, gene_c_tuple, attr_dict={
        RELATION: INCREASES,
    })

    graph.add_edge(gene_c_tuple, protein_g_tuple, attr_dict={
        RELATION: ASSOCIATION,
    })

    return graph


def make_graph_4():
    """
    A -> B
            B -| C
            B -| D
            B -| E
            B -| F
            B -> G

            B -> H
            B -| H

            B -> I
            B -- J

    """
    graph = BELGraph(
        name='PyBEL Tools Example Network 4',
        version='1.0.0',
        description='Example Network for PyBEL Tools Tests',
        authors='Daniel Domingo-Fernández and Charles Tapley Hoyt',
        contact='charles.hoyt@scai.fraunhofer.de',
    )

    graph.add_edge(protein_a_tuple, protein_b_tuple, attr_dict={
        RELATION: INCREASES,
    })

    graph.add_edge(protein_b_tuple, gene_c_tuple, attr_dict={
        RELATION: DECREASES,
    })

    graph.add_edge(protein_b_tuple, rna_d_tuple, attr_dict={
        RELATION: DECREASES,
    })

    graph.add_edge(protein_b_tuple, protein_e_tuple, attr_dict={
        RELATION: DECREASES,
    })

    graph.add_edge(protein_b_tuple, gene_f_tuple, attr_dict={
        RELATION: DECREASES,
    })

    graph.add_edge(protein_b_tuple, protein_g_tuple, attr_dict={
        RELATION: INCREASES,
    })

    graph.add_edge(protein_b_tuple, protein_h_tuple, attr_dict={
        RELATION: DECREASES,
    })

    graph.add_edge(protein_b_tuple, protein_h_tuple, attr_dict={
        RELATION: INCREASES,
    })

    graph.add_edge(protein_b_tuple, protein_i_tuple, attr_dict={
        RELATION: INCREASES,
    })

    graph.add_edge(protein_b_tuple, protein_j_tuple, attr_dict={
        RELATION: ASSOCIATION,
    })

    return graph


example_1 = make_graph_1()
example_2 = make_graph_2()
example_3 = make_graph_3()
example_4 = make_graph_4()


class ExampleNetworkMixin(unittest.TestCase):
    """A mixin that gives a class access to example networks"""

    def setUp(self):
        super(ExampleNetworkMixin, self).setUp()

        self.graph_1 = make_graph_1()
        self.network2 = make_graph_2()
        self.network3 = make_graph_3()
        self.network4 = make_graph_4()

def n():
    return str(uuid4())
