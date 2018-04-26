# -*- coding: utf-8 -*-

import sys

import click

from .algorithm import run_epicom
from ..neurommsig import get_ad_graph, get_ep_graph, get_pd_graph


@click.command()
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
def main(output):
    graphs = [
        get_ad_graph(),
        get_ep_graph(),
        get_pd_graph(),
    ]

    run_epicom(graphs, output)


if __name__ == '__main__':
    main()
