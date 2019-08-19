# -*- coding: utf-8 -*-

"""Utility for assigning pathways."""

import itertools as itt
import random
from collections import defaultdict

from pybel import BELGraph
from pybel.dsl import BaseAbundance, ListAbundance

__all__ = [
    'PathwayAssigner',
]


def is_hgnc(node):
    try:
        return node.namespace.lower() == 'hgnc'
    except AttributeError:
        return False


class PathwayAssigner:
    """A tool for assigning pathways to unannotated BEL edges."""

    def __init__(self, *, graph: BELGraph, managers):
        """Initialize the pathway assigner with several lookup dictionaries.

        :param managers: A ComPath manager or iterable of ComPath managers
        """
        self.graph = graph

        self.pathway_to_symbols = defaultdict(set)
        self.symbol_to_pathways = defaultdict(set)

        if not isinstance(managers, list):
            managers = []

        for manager in managers:
            self._add_manager(manager)

        # These won't be loaded more so convert to normal dicts
        self.pathway_to_symbols = dict(self.pathway_to_symbols)
        self.symbol_to_pathways = dict(self.symbol_to_pathways)

        self.pathway_to_key = defaultdict(set)
        self.key_to_pathway = defaultdict(set)

        self.pmid_to_pathway = defaultdict(set)
        self.pathway_to_pmid = defaultdict(set)

        self.double_annotated = defaultdict(lambda: defaultdict(list))

    def _add_manager(self, manager) -> None:
        """Add a ComPath manager."""
        for pathway in manager._query_pathway().all():
            db_id = getattr(pathway, f'{manager.module_name}_id')
            pathway_tuple = manager.module_name, db_id, pathway.name

            for protein in pathway.proteins:
                self.pathway_to_symbols[pathway_tuple].add(protein.hgnc_symbol)
                self.symbol_to_pathways[protein.hgnc_symbol].add(pathway_tuple)

    def to_file(self, tsv_path, rst_path):
        """Save results to files."""
        graph = self.graph

        with open(tsv_path, 'w') as file, open(rst_path, 'w') as log_file:
            print('database', 'pathway_id', 'pathway_name', 'key', 'bel', sep='\t', file=file)
            for (db, pathway_id, pathway), names_dict in self.double_annotated.items():
                title = f'{db}:{pathway_id} - {pathway}'
                print(title, file=log_file)
                print('=' * len(title), file=log_file)

                for node_key, keys_and_data in names_dict.items():
                    print('', file=log_file)
                    print(node_key, file=log_file)
                    print('-' * len(str(node_key)), file=log_file)
                    for u, v, key, data in keys_and_data:
                        print('-', key[:8], graph.edge_to_bel(u, v, data), file=log_file)
                        print(db, pathway_id, pathway, key, graph.edge_to_bel(u, v, data), sep='\t', file=file)

                print('', file=log_file)

    def summarize(self):
        """Print the summary of the annotations."""
        graph = self.graph

        annotated_edge_keys = set(itt.chain.from_iterable(self.pathway_to_key.values()))
        n_edges_annotated = len(annotated_edge_keys)

        print(f'{n_edges_annotated} ({n_edges_annotated / graph.number_of_edges():.2%}) '
              f'of {graph.number_of_edges()} edges were annotated')

        unannotated_edges = [
            (u, v, k, d)
            for u, v, k, d in graph.edges(data=True, keys=True)
            if k not in annotated_edge_keys
        ]

        print(f'There are {len(unannotated_edges)} unannotated edges')

        print('\nExamples of unannotated nodes:\n')
        for u, v, k, d in random.sample(unannotated_edges, 15):
            print(k[:8], graph.edge_to_bel(u, v, d))

        print()

        annotated_nodes = {
            node
            for u, v, k in graph.edges(keys=True)
            if k in annotated_edge_keys
            for node in (u, v)
        }

        n_nodes_annotated = len(annotated_nodes)

        print(f'{n_nodes_annotated} ({n_nodes_annotated / graph.number_of_nodes():.2%}) '
              f'of {graph.number_of_nodes()} nodes were annotated')

        unannotated_nodes = set(graph) - annotated_nodes
        print(f'There are {len(unannotated_nodes)} unannotated edges')

        print('\nExamples of unannotated nodes:\n')
        for node in random.sample(unannotated_nodes, 15):
            print(node)

    def annotate_gene_gene(self):
        """Annotate edges between gene or gene product nodes.

        1. `If` the subject and object in an edge are both in a canonical pathway, then the edge gets assigned to the
           pathway.
        2. `Else if` only one of the subject and the object in the edge have been assigned in the pathway:
          1. `If` the edge is an ontological edge, than add it to the pathway
          2. `If` there are other edges in the pathway mentioned in the same article, assign the edge to the pathway
          3. `Else` leave for manual curation
        3. `Else if` neither of the nodes are assigned to the pathway, but both nodes are connected to nodes in the
           pathway by directed edges, assign both edge to the pathway as well as incident edges
        4. `Else` the nodes don't get assigned to the pathway
        """
        graph = self.graph
        c = 0

        for u, v, k, d in graph.edges(keys=True, data=True):
            if not isinstance(u, BaseAbundance) or not isinstance(v, BaseAbundance):
                continue

            if not is_hgnc(u) or not is_hgnc(v):
                continue

            u_name, v_name = u.name, v.name

            for pathway_tuple, symbols in self.pathway_to_symbols.items():
                if u_name not in symbols or v_name not in symbols:
                    continue

                self.double_annotated[pathway_tuple][tuple(sorted([u_name, v_name]))].append((u, v, k, d))

                self.pathway_to_key[pathway_tuple].add(k)
                self.key_to_pathway[k].add(pathway_tuple)

                citation = d.get('citation')
                if citation is not None:
                    reference = citation['reference']
                    self.pmid_to_pathway[reference].add(pathway_tuple)
                    self.pathway_to_pmid[pathway_tuple].add(reference)

                c += 1

        return c

    def annotate_gene_other(self):
        """Annotate edges between gene or gene product nodes and chemicals / biological processes / diseases.

        If an entity is related to a gene in a pathway, then that edge gets annotated to the pathway
        """
        graph = self.graph
        c = 0

        for u, v, k, d in graph.edges(keys=True, data=True):
            if not isinstance(u, BaseAbundance) or not isinstance(v, BaseAbundance):
                continue

            if is_hgnc(u) and not is_hgnc(v):
                gene_name = u.name
                other_name = v.name
            elif not is_hgnc(u) and is_hgnc(v):
                gene_name = v.name
                other_name = u.name
            else:
                continue

            for pathway_tuple, symbols in self.pathway_to_symbols.items():
                if gene_name not in symbols:
                    continue

                self.double_annotated[pathway_tuple][tuple(sorted([gene_name, other_name]))].append((u, v, k, d))

                self.pathway_to_key[pathway_tuple].add(k)
                self.key_to_pathway[k].add(pathway_tuple)

                citation = d.get('citation')
                if citation is not None:
                    reference = citation['reference']
                    self.pmid_to_pathway[reference].add(pathway_tuple)
                    self.pathway_to_pmid[pathway_tuple].add(reference)

                c += 1

        return c

    def annotate_by_document(self):
        """Annotate edges with a gene or gene product nodes that has other annotated edges in its original document.

        If an edge has only one node that appears in a pathway, but that pathway has already been mentioned in the
        paper, then it gets annotated to that pathway too.
        """
        graph = self.graph
        c = 0

        for u, v, k, d in graph.edges(keys=True, data=True):
            citation = d.get('citation')
            if citation is None:
                continue

            reference = citation['reference']
            pathway_tuples = self.pmid_to_pathway[reference]

            if is_hgnc(u) and not is_hgnc(v):
                gene_name = u.name
            elif not is_hgnc(u) and is_hgnc(v):
                gene_name = v.name
            else:
                continue

            if gene_name not in self.symbol_to_pathways:
                continue

            for pathway_tuple in pathway_tuples:
                if pathway_tuple not in self.symbol_to_pathways[gene_name]:
                    continue

                self.double_annotated[pathway_tuple][gene_name].append((u, v, k, d))

                self.pathway_to_key[pathway_tuple].add(k)
                self.key_to_pathway[k].add(pathway_tuple)

                c += 1

        return c

    def annotate_complexes(self):
        """Annotated complex nodes.

        If two or more members of a complex are in a pathway, then the whole complex and all of its partOf
        relationships will get assigned to that pathway.
        """
        graph = self.graph
        c = 0

        for node in graph:
            if not isinstance(node, ListAbundance):
                continue

            hgnc_count = sum(
                member.namespace.lower() == 'hgnc'
                for member in node.members
                if isinstance(member, BaseAbundance)
            )

            if 0 == hgnc_count:
                continue

            for pathway_tuple, symbols in self.pathway_to_symbols.items():
                in_count = sum(
                    member.name in symbols
                    for member in node.members
                    if isinstance(member, BaseAbundance) and member.namespace.lower() == 'hgnc'
                )

                do_it = (
                    (1 == hgnc_count and 1 == in_count)  # Other stuff going on, lets do it
                    or 2 <= in_count  # enough is going on

                )

                if not do_it:
                    continue

                for u, v, k, d in graph.edges(node, keys=True, data=True):
                    self.double_annotated[pathway_tuple][node].append((u, v, k, d))

                    self.pathway_to_key[pathway_tuple].add(k)
                    self.key_to_pathway[k].add(pathway_tuple)

                    c += 1

        return c

    # TODO add FamPlex hierarchy resolution
    # TODO add orthology resolution
    # TODO add partOf relationship resolution
