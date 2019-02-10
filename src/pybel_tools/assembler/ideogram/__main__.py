# -*- coding: utf-8 -*-

import os

import click

from .assembler import HERE, to_html


@click.command()
@click.option('-o', '--output', type=click.File('w'), default=os.path.join(HERE, 'test.html'))
def main(output):
    """Output the HBP knowledge graph to the desktop"""
    from hbp_knowledge import get_graph
    graph = get_graph()
    text = to_html(graph)
    print(text, file=output)


if __name__ == '__main__':
    main()
