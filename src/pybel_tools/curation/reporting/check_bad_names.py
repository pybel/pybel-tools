# -*- coding: utf-8 -*-

"""This script checks which of the BEL graphs have unhelpful names."""

import re

import click
from hbp_knowledge import get_graphs

BAD_NAME = re.compile(r'.*\d{4}$')


@click.command()
def main():
    """Check which graphs' names end with a year."""
    for name, graph in get_graphs().items():
        if BAD_NAME.match(graph.name):
            click.echo(f'{graph.name} {graph.path}')


if __name__ == '__main__':
    main()
