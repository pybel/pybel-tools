# -*- coding: utf-8 -*-

"""Utility for assigning pathways."""

import itertools as itt
import json
import os
import random
from collections import defaultdict
from typing import List, Optional

import bio2bel_hgnc
import bio2bel_mgi
import bio2bel_rgd
from pybel import BELGraph
from pybel.dsl import BaseAbundance, ListAbundance

__all__ = [
    'PathwayAssigner',
]


class PathwayAssigner:
    """A tool for assigning pathways to unannotated BEL edges."""

    def __init__(
        self,
        *,
        graph: BELGraph,
        managers: List,
        mgi_cache_path: Optional[str] = None,
        rgd_cache_path: Optional[str] = None,
    ):
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

        # Prepare RGD mappings
        if rgd_cache_path is not None and os.path.exists(rgd_cache_path):
            with open(rgd_cache_path) as file:
                self.rgd_symbol_to_hgnc_symbol = json.load(file)
        else:
            self.rgd_symbol_to_hgnc_symbol = {}
            rgd_symbol_to_id = {}
            rgd_id_to_symbol = {}

            rgd_manager = bio2bel_rgd.Manager()
            rat_genes = rgd_manager.list_genes()
            for rat_gene in rat_genes:
                rgd_symbol_to_id[rat_gene.symbol] = rat_gene.rgd_id
                rgd_id_to_symbol[rat_gene.rgd_id] = rat_gene.symbol

            hgnc_manager = bio2bel_hgnc.Manager()
            genes = hgnc_manager.list_human_genes()
            for gene in genes:
                for rgd_id in gene.rgds:
                    rgd_id = str(rgd_id)
                    rgd_symbol = rgd_id_to_symbol.get(rgd_id)
                    if rgd_symbol is None:
                        print(f'could not find rgd:{rgd_id}')
                        continue
                    self.rgd_symbol_to_hgnc_symbol[rgd_symbol] = gene.symbol

            if rgd_cache_path is not None:
                with open(rgd_cache_path, 'w') as file:
                    json.dump(self.rgd_symbol_to_hgnc_symbol, file, indent=2)

        # Prepare MGI mappings
        if mgi_cache_path is not None and os.path.exists(mgi_cache_path):
            with open(mgi_cache_path) as file:
                self.mgi_symbol_to_hgnc_symbol = json.load(file)
        else:
            self.mgi_symbol_to_hgnc_symbol = {}
            mgi_symbol_to_id = {}
            mgi_id_to_symbol = {}
            mgi_to_hgnc = {}

            mgi_manager = bio2bel_mgi.Manager()
            mouse_genes = mgi_manager.get_genes()
            for mouse_gene in mouse_genes:
                mgi_symbol_to_id[mouse_gene.symbol] = mouse_gene.mgi_id
                mgi_id_to_symbol[mouse_gene.mgi_id] = mouse_gene.symbol

            hgnc_manager = bio2bel_hgnc.Manager()
            genes = hgnc_manager.list_human_genes()
            for gene in genes:
                for mgi_id in gene.mgds:
                    mgi_id = f'MGI:{mgi_id}'
                    mgi_to_hgnc[mgi_id] = (gene.identifier, gene.symbol)

                    mgi_symbol = mgi_id_to_symbol.get(mgi_id)
                    if mgi_symbol is None:
                        print(f'could not find {mgi_id}')
                        continue
                    self.mgi_symbol_to_hgnc_symbol[mgi_symbol] = gene.symbol

            if mgi_cache_path is not None:
                with open(mgi_cache_path, 'w') as file:
                    json.dump(self.mgi_symbol_to_hgnc_symbol, file, indent=2)

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
        print(f'There are {len(unannotated_nodes)} unannotated nodes')

        print('\nExamples of unannotated nodes:\n')
        for node in random.sample(unannotated_nodes, 15):
            print(node)

    def get_gene(self, node: BaseAbundance) -> Optional[str]:
        """Get or map the name to HGNC gene symbol of the node, if possible."""
        try:
            namespace = node.namespace.lower()
        except AttributeError:
            return

        if namespace == 'hgnc':
            return node.name

        if namespace == 'rgd':
            return self.rgd_symbol_to_hgnc_symbol.get(node.name)

        if namespace == 'mgi':
            return self.mgi_symbol_to_hgnc_symbol.get(node.name)

    def annotate_gene_gene(self):
        """Annotate edges between gene or gene product nodes.

        1. Identify if subject and object are both gene nodes. If they are orthologs, try and map them to HGNC.
        2. `If` the subject and object in an edge are both in a canonical pathway, then the edge gets assigned to the
           pathway.
        3. `Else if` only one of the subject and the object in the edge have been assigned in the pathway:
          1. `If` the edge is an ontological edge, than add it to the pathway
          2. `If` there are other edges in the pathway mentioned in the same article, assign the edge to the pathway
          3. `Else` leave for manual curation
        4. `Else if` neither of the nodes are assigned to the pathway, but both nodes are connected to nodes in the
           pathway by directed edges, assign both edge to the pathway as well as incident edges
        5. `Else` the nodes don't get assigned to the pathway
        """
        graph = self.graph
        c = 0

        for u, v, k, d in graph.edges(keys=True, data=True):
            if not isinstance(u, BaseAbundance) or not isinstance(v, BaseAbundance):
                continue

            u_name, v_name = self.get_gene(u), self.get_gene(v)
            if u_name is None or v_name is None:
                continue

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

        1. Identify if subject or object are a gene nodes. If they are orthologs, try and map them to HGNC.
        2. If an entity is related to a gene in a pathway, then that edge gets annotated to the pathway
        """
        graph = self.graph
        c = 0

        for u, v, k, d in graph.edges(keys=True, data=True):
            if not isinstance(u, BaseAbundance) or not isinstance(v, BaseAbundance):
                continue

            u_name, v_name = self.get_gene(u), self.get_gene(v)
            if u_name and v_name is None:
                gene_name = u_name
                other_name = v.name
            elif u_name is None and v_name:
                gene_name = v_name
                other_name = u.name
            else:
                continue

            try:
                ordering = tuple(sorted([gene_name, other_name]))
            except TypeError:
                print('Gene', gene_name)
                print('Other', other_name)
                continue

            for pathway_tuple, symbols in self.pathway_to_symbols.items():
                if gene_name not in symbols:
                    continue

                self.double_annotated[pathway_tuple][ordering].append((u, v, k, d))

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

            u_name, v_name = self.get_gene(u), self.get_gene(v)
            if u_name and v_name is None:
                gene_name = u_name
            elif u_name is None and v_name:
                gene_name = v_name
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

            mapped_names = []
            for member in node.members:
                if not isinstance(member, BaseAbundance):
                    continue
                name = self.get_gene(member)
                if name is not None:
                    mapped_names.append(name)

            if not mapped_names:
                continue

            for pathway_tuple, symbols in self.pathway_to_symbols.items():
                in_count = sum(
                    name in symbols
                    for name in mapped_names
                )
                should_annotate_complex = (
                    (1 == len(mapped_names) and 1 == in_count)  # Other stuff going on, let's do it
                    or 2 <= in_count  # enough is going on, let's do it
                )
                if not should_annotate_complex:
                    continue
                # do it
                # do it
                # do it
                for u, v, k, d in graph.edges(node, keys=True, data=True):
                    self.double_annotated[pathway_tuple][node].append((u, v, k, d))

                    self.pathway_to_key[pathway_tuple].add(k)
                    self.key_to_pathway[k].add(pathway_tuple)

                    c += 1

        return c

    # TODO add FamPlex hierarchy resolution
    # TODO add partOf relationship resolution
