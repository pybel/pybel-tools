# -*- coding: utf-8 -*-

"""Tools for identifying common HGVS mistakes in BEL graphs."""

from typing import Iterable, Mapping, Optional, TextIO, Tuple

from pybel import BELGraph
from pybel.dsl import CentralDogma, Hgvs
from pybel.language import amino_acid_dict, dna_nucleotide_labels
from pybel.struct import has_hgvs

__all__ = [
    'summarize_graphs_variants',
    'iterate_variants',
]


def summarize_graphs_variants(graphs: Mapping[str, BELGraph], file: Optional[TextIO] = None) -> None:
    """Summarize a group of BEL graphs."""
    for graph_name, graph in graphs.items():
        hgvss = list(iterate_variants(graph))
        if not hgvss:
            continue
        print(graph_name, graph.authors, file=file)
        for node, old, new in hgvss:
            if new:
                print(f'  {node.name}: {old} -> {new}', file=file)
            else:
                print(f'  {node.name}: {old}', file=file)


def iterate_variants(graph: BELGraph) -> Iterable[Tuple[CentralDogma, str, Optional[str]]]:
    """Yield nodes with HGVS variants that might need to be fixed."""
    for node in graph:
        if not (isinstance(node, CentralDogma) and has_hgvs(node)):
            continue
        for variant in node.variants:
            if not isinstance(variant, Hgvs):
                continue
            hgvs_string = variant['identifier']
            if is_common_hgvs(hgvs_string):
                continue

            if is_one_letter_protein_substitution(hgvs_string):
                yield node, hgvs_string, upgrade_one_letter_protein_substitution(hgvs_string)
            elif is_one_letter_protein_truncation(hgvs_string):
                yield node, hgvs_string, upgrade_one_letter_protein_truncation(hgvs_string)
            else:
                yield node, hgvs_string, None


def _find_replace(path: str, old: str, new: str) -> None:
    """Find and replace a given string in a file with."""
    with open(path) as file:
        lines = [
            line.strip().replace(old, new)
            for line in file
        ]
    with open(path, 'w') as file:
        for line in lines:
            print(line, file=file)


def is_common_hgvs(hgvs_string: str) -> bool:
    """Check if the HGVS string is a common case."""
    return (
        is_protein_substitution(hgvs_string)
        or hgvs_string == '?'
        or is_protein_deletion(hgvs_string)
        or is_protein_truncation(hgvs_string)
        or is_gene_substitution(hgvs_string)
    )


def is_protein_substitution(hgvs_string: str) -> bool:
    """Check if the HGVS string is a protein substitution."""
    return (
        hgvs_string.startswith('p.')
        and hgvs_string[2:5] in amino_acid_dict.values()
        and hgvs_string[-3:] in amino_acid_dict.values()
    )


def is_protein_deletion(hgvs_string: str) -> bool:
    """Check if the HGVS string is a protein deletion."""
    return (
        hgvs_string.startswith('p.')
        and hgvs_string[2:5] in amino_acid_dict.values()
        and hgvs_string[-3:] == 'del'
    )


def is_protein_truncation(hgvs_string: str) -> bool:
    """Check if the HGVS string is a protein truncation."""
    return (
        hgvs_string.startswith('p.')
        and hgvs_string[2:5] in amino_acid_dict.values()
        and hgvs_string[-1:] == '*'
    )


def is_one_letter_protein_substitution(hgvs_string: str) -> bool:
    """Check if the HGVS string is a malformed protein substitution."""
    return (
        hgvs_string.startswith('p.')
        and hgvs_string[2] in amino_acid_dict
        and hgvs_string[-1] in amino_acid_dict
    )


def is_one_letter_protein_truncation(hgvs_string: str) -> bool:
    """Check if the HGVS string is a malformed protein truncation."""
    return (
        hgvs_string.startswith('p.')
        and hgvs_string[2] in amino_acid_dict
        and hgvs_string[-1] == '*'
    )


def is_gene_substitution(hgvs_string: str) -> bool:
    """Check if the HGVS string is a genetic substitution."""
    return (
        hgvs_string.startswith('c.')
        and hgvs_string[-3] in dna_nucleotide_labels
        and hgvs_string[-2] == '>'
        and hgvs_string[-1] in dna_nucleotide_labels
    )


def upgrade_one_letter_protein_substitution(hgvs_string: str) -> str:
    """Upgrade a malformed one letter protein substitution."""
    return f'p.{amino_acid_dict[hgvs_string[2]]}{hgvs_string[3:-1]}{amino_acid_dict[hgvs_string[-1]]}'


def upgrade_one_letter_protein_truncation(hgvs_string: str) -> str:
    """Upgrade a malformed one letter protein truncation."""
    return f'p.{amino_acid_dict[hgvs_string[2]]}{hgvs_string[3:-1]}*'


def _main():
    from hbp_knowledge import get_graphs
    summarize_graphs_variants(get_graphs())


if __name__ == '__main__':
    _main()
