# -*- coding: utf-8 -*-

import os
import tempfile
import unittest

from pybel import BELGraph
from pybel.dsl import gene, protein, rna
from pybel.manager import Manager
from pybel.testing.utils import n

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


def make_graph_1() -> BELGraph:
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

    graph.add_increases(
        protein_a,
        protein_b,
        citation='1',
        evidence='Evidence 1',
        annotations={'Annotation': 'foo'}
    )

    graph.add_increases(
        rna_d,
        protein_a,
        citation='2',
        evidence='Evidence 2',
        annotations={'Annotation': 'foo'}
    )

    graph.add_decreases(
        gene_c,
        protein_b,
        citation='3',
        evidence='Evidence 3',
        annotations={'Annotation': 'foo'}
    )

    return graph


def make_graph_2() -> BELGraph:
    """Make an example graph."""
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

    graph.add_increases(
        protein_e,
        protein_b,
        citation='1',
        evidence='Evidence 1',
        annotations={'Annotation': 'foo'},
    )

    graph.add_increases(
        gene_f,
        protein_e,
        citation='2',
        evidence='Evidence 2',
        annotations={'Annotation': 'foo2'}
    )

    return graph


def make_graph_3() -> BELGraph:
    """Make an example graph.

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

    graph.add_increases(protein_a, protein_b, n(), n())
    graph.add_decreases(protein_b, gene_c, n(), n())
    graph.add_decreases(rna_d, gene_f, n(), n())
    graph.add_increases(protein_e, gene_f, n(), n())
    graph.add_increases(gene_f, gene_c, n(), n())
    graph.add_association(gene_c, protein_g, n(), n())

    return graph


def make_graph_4() -> BELGraph:
    """Make an example graph.

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

    graph.add_increases(protein_a, protein_b, n(), n())
    graph.add_decreases(protein_b, gene_c, n(), n())
    graph.add_decreases(protein_b, rna_d, n(), n())
    graph.add_decreases(protein_b, protein_e, n(), n())
    graph.add_decreases(protein_b, gene_f, n(), n())
    graph.add_increases(protein_b, protein_g, n(), n())
    graph.add_decreases(protein_b, protein_h, n(), n())
    graph.add_increases(protein_b, protein_h, n(), n())
    graph.add_increases(protein_b, protein_i, n(), n())
    graph.add_association(protein_b, protein_j, n(), n())

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
