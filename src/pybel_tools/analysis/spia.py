# -*- coding: utf-8 -*-

"""An exporter for signaling pathway impact analysis (SPIA) described by [Tarca2009]_.

.. [Tarca2009] Tarca, A. L., *et al* (2009). `A novel signaling pathway impact analysis
               <https://doi.org/10.1093/bioinformatics/btn577>`_. Bioinformatics, 25(1), 75–82.

To run this module on an arbitrary BEL graph, use the command ``python -m pybel_tools.analysis.spia``.

.. seealso:: https://bioconductor.org/packages/release/bioc/html/SPIA.html
"""

from typing import Dict, Set, Mapping

import click
import pandas as pd
from itertools import product

from pybel import BELGraph
from pybel.cli import graph_pickle_argument
from pybel.constants import ASSOCIATION, CAUSAL_DECREASE_RELATIONS, CAUSAL_INCREASE_RELATIONS, IDENTIFIER, NAME
from pybel.dsl import CentralDogma, Gene, ListAbundance, ProteinModification, Rna

__all__ = [
    'run',
    'bel_to_spia_matrices',
    'spia_matrices_to_excel',
]

KEGG_RELATIONS = {
    "activation",
    "compound",
    "binding_association",
    "expression",
    "inhibition",
    "activation_phosphorylation",
    "phosphorylation",
    "inhibition_phosphorylation",
    "inhibition_dephosphorylation",
    "dissociation",
    "dephosphorylation",
    "activation_dephosphorylation",
    "state change",
    "activation_indirect effect",
    "inhibition_ubiquination",
    "ubiquination",
    "expression_indirect effect",
    "inhibition_indirect effect",
    "repression",
    "dissociation_phosphorylation",
    "indirect effect_phosphorylation",
    "activation_binding_association",
    "indirect effect",
    "activation_compound",
    "activation_ubiquination"
}


def run(graph: BELGraph, path: str) -> None:
    """Run the SPIA pipeline and export to an excel at the given location."""
    spia_matrices = bel_to_spia_matrices(graph)
    spia_matrices_to_excel(spia_matrices, path)


def bel_to_spia_matrices(graph: BELGraph) -> Mapping[str, pd.DataFrame]:
    """Create an excel sheet ready to be used in SPIA software.

    :param graph: BELGraph
    :return: dictionary with matrices
    """
    index_nodes = get_matrix_index(graph)
    spia_matrices = build_spia_matrices(index_nodes)

    for sub, obj, data in graph.edges(data=True):
        # Both nodes are CentralDogma abundances
        if isinstance(sub, CentralDogma) and isinstance(obj, CentralDogma):
            # Update matrix dict
            update_spia_matrices(spia_matrices, sub, obj, data)

        # Subject is CentralDogmaAbundance and node is ListAbundance
        elif isinstance(sub, CentralDogma) and isinstance(obj, ListAbundance):
            # Add a relationship from subject to each of the members in the object
            for node in obj.members:
                # Skip if the member is not in CentralDogma
                if not isinstance(node, CentralDogma):
                    continue

                update_spia_matrices(spia_matrices, sub, node, data)

        # Subject is ListAbundance and node is CentralDogmaAbundance
        elif isinstance(sub, ListAbundance) and isinstance(obj, CentralDogma):
            # Add a relationship from each of the members of the subject to the object
            for node in sub.members:
                # Skip if the member is not in CentralDogma
                if not isinstance(node, CentralDogma):
                    continue

                update_spia_matrices(spia_matrices, node, obj, data)

        # Both nodes are ListAbundance
        elif isinstance(sub, ListAbundance) and isinstance(obj, ListAbundance):
            for sub_member, obj_member in product(sub.members, obj.members):
                # Update matrix if both are CentralDogma
                if isinstance(sub_member, CentralDogma) and isinstance(obj_member, CentralDogma):
                    update_spia_matrices(spia_matrices, sub_member, obj_member, data)

        # else Not valid edge
    return spia_matrices


def get_matrix_index(graph: BELGraph) -> Set[str]:
    """Return set of HGNC names from Proteins/Rnas/Genes/miRNA, nodes that can be used by SPIA."""
    # TODO: Using HGNC Symbols for now
    return {
        node.name
        for node in graph
        if isinstance(node, CentralDogma) and node.namespace.upper() == 'HGNC'
    }


def build_spia_matrices(nodes: Set[str]) -> Dict[str, pd.DataFrame]:
    """Build an adjacency matrix for each KEGG relationship and return in a dictionary.

    :param nodes: A set of HGNC gene symbols
    :return: Dictionary of adjacency matrix for each relationship
    """
    nodes = list(sorted(nodes))
    return {
        relation: pd.DataFrame(0, index=nodes, columns=nodes)
        for relation in KEGG_RELATIONS
    }


def update_spia_matrices(spia_matrices: Dict[str, pd.DataFrame], sub: CentralDogma, obj: CentralDogma, data: Dict):
    """Populate the adjacency matrix.

    :param spia_matrices:
    :param sub: bel subject
    :param obj: bel object
    :param data: edge data
    """
    if sub.namespace.upper() != 'HGNC' and obj.namespace.upper() != 'HGNC':
        return

    subject_name = sub.name
    object_name = obj.name
    relation = data['relation']

    if relation in CAUSAL_INCREASE_RELATIONS:
        # If it has pmod check which one and add it to the corresponding matrix
        if obj.variants and any(isinstance(variant, ProteinModification) for variant in obj.variants):
            for variant in obj.variants:
                if variant[IDENTIFIER][NAME] == "Ub":
                    spia_matrices["activation_ubiquination"][subject_name][object_name] = 1

                elif variant[IDENTIFIER][NAME] == "Ph":
                    spia_matrices["activation_phosphorylation"][subject_name][object_name] = 1

        # Normal increase, add activation
        elif isinstance(obj, (Gene, Rna)):
            spia_matrices['expression'][subject_name][object_name] = 1

        else:
            spia_matrices['activation'][subject_name][object_name] = 1

    elif relation in CAUSAL_DECREASE_RELATIONS:
        # If it has pmod check which one and add it to the corresponding matrix
        if obj.variants and any(isinstance(variant, ProteinModification) for variant in obj.variants):
            for variant in obj.variants:
                if variant[IDENTIFIER][NAME] == "Ub":
                    spia_matrices['inhibition_ubiquination'][subject_name][object_name] = 1

                elif variant[IDENTIFIER][NAME] == "Ph":
                    spia_matrices["inhibition_phosphorylation"][subject_name][object_name] = 1

        # Normal decrease, check which matrix
        elif isinstance(obj, (Gene, Rna)):
            spia_matrices["repression"][subject_name][object_name] = 1

        else:
            spia_matrices["inhibition"][subject_name][object_name] = 1

    elif relation == ASSOCIATION:
        spia_matrices["binding_association"][subject_name][object_name] = 1


def spia_matrices_to_excel(spia_matrices: Mapping[str, pd.DataFrame], path: str) -> None:
    """Export a SPIA data dictionary into an Excel sheet at the given path.

    .. note::

        # The R import should add the values:
        # ["nodes"] from the columns
        # ["title"] from the name of the file
        # ["NumberOfReactions"] set to "0"
    """
    writer = pd.ExcelWriter('{}.xlsx'.format(path), engine='xlsxwriter')

    for relation, df in spia_matrices.items():
        df.to_excel(writer, sheet_name=relation)

    # Save excel
    writer.save()


@click.command()
@graph_pickle_argument
@click.argument('output')
def main(graph: BELGraph, output: str):
    """Export the graph to a SPIA Excel sheet."""
    run(graph, output)


if __name__ == '__main__':
    main()
