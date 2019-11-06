# -*- coding: utf-8 -*-

"""Summarize how much effort each curator has put."""

from collections import Counter, defaultdict
from typing import Mapping

import click

from pybel import BELGraph


def summarize_effort_per_graph(graphs: Mapping[str, BELGraph]):
    """Summarize the effort per curator."""
    r = defaultdict(list)
    for name, graph in graphs.items():
        if ',' in graph.authors:
            first_author = graph.authors.split(',')[0].strip()
        elif ' and ' in graph.authors:
            i = graph.authors.find(' and ')
            first_author = graph.authors[:i]
        else:  # sole author
            first_author = graph.authors

        r[first_author].append(graph)

    max_author_width = max(map(len, r))

    click.echo('= Author Graphs =')
    graph_counter = Counter({
        author: len(graphs)
        for author, graphs in r.items()
    })
    for author, graph_count in graph_counter.most_common():
        click.echo(f'{author:{max_author_width}} {graph_count:>}')

    click.echo('\n= Author Edges =')
    edge_counter = Counter({
        author: sum(graph.number_of_edges() for graph in graphs)
        for author, graphs in r.items()
    })
    for author, edge_count in edge_counter.most_common():
        click.echo(f'{author:{max_author_width}} {edge_count:>}')


@click.command()
def main():
    """Summarize the effort per curator."""
    import hbp_knowledge

    graphs: Mapping[str, BELGraph] = hbp_knowledge.get_graphs()
    summarize_effort_per_graph(graphs)


if __name__ == '__main__':
    main()
