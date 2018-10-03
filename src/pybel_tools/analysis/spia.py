# -*- coding: utf-8 -*-

"""This module contains the BEL exporter to SPIA format.
(see: https://bioconductor.org/packages/release/bioc/html/SPIA.html)
"""
import typing
from itertools import product

import pandas as pd

from pybel import BELGraph
from pybel.constants import ASSOCIATION, INCREASES, DECREASES, DIRECTLY_DECREASES, DIRECTLY_INCREASES, IDENTIFIER, NAME
from pybel.dsl.node_classes import BaseAbundance, CentralDogma, Gene, ListAbundance, Rna, ProteinModification
from pybel_tools.constants import KEGG_RELATIONS


def get_matrix_index(graph: BELGraph) -> typing.Set:
    """Return set of HGNC names from Proteins/Rnas/Genes/miRNA, nodes that can be used by SPIA."""
    # TODO: Using HGNC Symbols for now
    return {
        node.name
        for node in graph.nodes()
        if isinstance(node, BaseAbundance) and node.namespace == 'HGNC'
    }


def build_matrices(nodes: typing.Set) -> typing.Dict[str, pd.DataFrame]:
    """Build adjacency matrices for each KEGG relationships.

    :param nodes: nodes to be included in the matrix
    :return: Dictionary of adjency matrix for each relationship
    """
    return {
        relation: pd.DataFrame(0, index=nodes, columns=nodes)
        for relation in KEGG_RELATIONS
    }


def update_matrix(matrix_dict, sub, obj, data):
    """Populate adjacency matrix.

    :param Dict[str, pd.DataFrame] matrix_dict:
    :param sub: bel subject
    :param obj: bel object
    :param dict data: edge data
    """
    if sub.namespace == 'HGNC' or obj.namespace == 'HGNC':
        subject_name = sub.name
        object_name = obj.name

        relation = data['relation']

        if relation in {DIRECTLY_INCREASES, INCREASES}:

            # If it has pmod check which one and add it to the corresponding matrix
            if obj.variants and any(isinstance(variant, ProteinModification) for variant in obj.variants):

                for variant in obj.variants:

                    if variant[IDENTIFIER][NAME] == "Ub":
                        matrix_dict["activation_ubiquination"][subject_name][object_name] = 1

                    elif variant[IDENTIFIER][NAME] == "Ph":
                        matrix_dict["activation_phosphorylation"][subject_name][object_name] = 1

            # Normal increase, add activation
            else:

                if isinstance(obj, (Gene, Rna)):
                    matrix_dict['expression'][subject_name][object_name] = 1

                else:
                    matrix_dict['activation'][subject_name][object_name] = 1

        if relation in {DIRECTLY_DECREASES, DECREASES}:

            # If it has pmod check which one and add it to the corresponding matrix
            if obj.variants and any(isinstance(variant, ProteinModification) for variant in obj.variants):

                for variant in obj.variants:

                    if variant[IDENTIFIER][NAME] == "Ub":
                        matrix_dict['inhibition_ubiquination'][subject_name][object_name] = 1

                    elif variant[IDENTIFIER][NAME] == "Ph":
                        matrix_dict["inhibition_phosphorylation"][subject_name][object_name] = 1

            # Normal decrease, check which matrix
            else:

                if isinstance(obj, Gene) or isinstance(obj, Rna):
                    matrix_dict["repression"][subject_name][object_name] = 1

                else:
                    matrix_dict["inhibition"][subject_name][object_name] = 1

        if relation == ASSOCIATION:
            matrix_dict["binding_association"][subject_name][object_name] = 1


def bel_to_spia(graph: BELGraph) -> typing.Dict:
    """Create excel sheet ready to be used in SPIA software.

    :param graph: BELGraph
    :return: dictionary with matrices
    """
    index_nodes = get_matrix_index(graph)
    matrix_dict = build_matrices(index_nodes)

    for sub, obj, data in graph.edges(data=True):

        # Both nodes are CentralDogma abundances
        if issubclass(type(sub), CentralDogma) and issubclass(type(obj), CentralDogma):
            # Update matrix dict
            update_matrix(matrix_dict, sub, obj, data)

        # Subject is CentralDogmaAbundance and node is ListAbundance
        elif issubclass(type(sub), CentralDogma) and issubclass(type(obj), ListAbundance):
            # Add a relationship from subject to each of the members in the object
            for node in obj.members:

                # Skip if the member is not in CentralDogma
                if not issubclass(type(node), CentralDogma):
                    continue

                update_matrix(matrix_dict, sub, node, data)

        # Subject is ListAbundance and node is CentralDogmaAbundance
        elif issubclass(type(sub), ListAbundance) and issubclass(type(obj), CentralDogma):
            # Add a relationship from each of the members of the subject to the object
            for node in sub.members:

                # Skip if the member is not in CentralDogma
                if not issubclass(type(node), CentralDogma):
                    continue

                update_matrix(matrix_dict, node, obj, data)

        # Both nodes are ListAbundance
        elif issubclass(type(sub), ListAbundance) and issubclass(type(obj), ListAbundance):
            for sub_member, obj_member in product(sub.members, obj.members):

                # Update matrix if both are CentralDogma
                if issubclass(type(sub_member), CentralDogma) and issubclass(type(obj_member), CentralDogma):
                    update_matrix(matrix_dict, sub_member, obj_member, data)

        # else Not valid edge
    return matrix_dict


def spia_to_excel(spia_data_dict: typing.Dict[str, pd.DataFrame], file_name: str):
    """Export SPIA data dictionary into an excel sheet.

    :param spia_data_dict: data coming from bel_to_spia
    :param file_name: file name
    """
    writer = pd.ExcelWriter('{}.xlsx'.format(file_name), engine='xlsxwriter')

    for relation, df in spia_data_dict.items():
        df.to_excel(writer, sheet_name=relation)

    # The R import should add the values:
    # ["nodes"] from the columns
    # ["title"] from the name of the file
    # ["NumberOfReactions"] set to "0"

    # Save excel
    writer.save()
