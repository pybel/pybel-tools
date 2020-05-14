# -*- coding: utf-8 -*-

"""Build a network of citations connected by co-occurrence of entities."""

from collections import Counter, defaultdict

import click
import itertools as itt
import networkx as nx
from tqdm import tqdm
from typing import TextIO

from pybel import BELGraph, Manager
from pybel.cli import connection_option, graph_pickle_argument
from pybel.constants import CITATION, CITATION_DB, CITATION_DB_NAME, CITATION_IDENTIFIER, CITATION_TYPE_PUBMED
from pybel.manager.citation_utils import enrich_pubmed_citations


@click.command()
@connection_option
@graph_pickle_argument
@click.option('-o', '--output', type=click.File('w'), required=True)
@click.option('-t', '--threshold', type=int, default=1)
def main(connection: str, graph: BELGraph, output: TextIO, threshold):
    """Build a citation network from the graph."""
    enrich_pubmed_citations(Manager(connection=connection), graph)
    citation_network = make_citation_network(graph, threshold=threshold)
    print('Source', 'Source Title', 'Target', 'Target Title', 'Shared', sep='\t', file=output)
    for u, v, d in citation_network.edges(data=True):
        print(
            u,
            citation_network.nodes[u]['title'],
            v,
            citation_network.nodes[v]['title'],
            d['weight'],
            sep='\t',
            file=output,
        )


def make_citation_network(bel_graph: BELGraph, threshold: int = 0) -> nx.Graph:
    """Make a citation network from the BEL graph based on which statements occur in multiple sourves."""
    dd = defaultdict(set)
    names = {}
    for u, v, k, d in bel_graph.edges(keys=True, data=True):
        citation = d.get(CITATION)
        if citation is None or citation[CITATION_DB] != CITATION_TYPE_PUBMED:
            continue
        reference = citation[CITATION_IDENTIFIER]
        dd[reference].update((u, v))
        names[reference] = citation.get(CITATION_DB_NAME)

    all_nodes = set(itt.chain.from_iterable(dd.values()))

    iterator = itt.product(all_nodes, itt.combinations(dd.items(), r=2))
    iterator = tqdm(iterator, total=len(all_nodes) * (len(dd) ** 2))
    c = Counter(
        (c1, c2)
        for node, ((c1, c1_values), (c2, c2_values)) in iterator
        if node in c1_values and node in c2_values
    )

    rv = nx.Graph()
    for (c1, c2), weight in c.items():
        if weight >= threshold:
            rv.add_edge(c1, c2, weight=weight)

    for reference, title in names.items():
        rv.nodes[reference]['title'] = title

    return rv


if __name__ == '__main__':
    main()
