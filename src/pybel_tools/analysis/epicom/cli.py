# -*- coding: utf-8 -*-

import os
import sys

import click

from pybel.cli import graph_pickle_argument
from .algorithm import multi_run_epicom, run_epicom
from ..neurommsig import get_ad_graph, get_ep_graph, get_pd_graph


directory_option = click.option('-d', '--directory', default=os.getcwd())

@click.group()
def main():
    """Run EpiCom Reloaded."""


@main.command()
@graph_pickle_argument
@directory_option
def run(graph, directory):
    """Run on an arbitrary graph."""
    run_epicom(graph, directory)


@main.command()
@directory_option
def ad(directory):
    """Run on the AD graph."""
    graph = get_ad_graph()
    run_epicom(graph, directory)


@main.command()
@directory_option
def pd(directory):
    """Run on the PD graph."""
    graph = get_pd_graph()
    run_epicom(graph, directory)


@main.command()
@directory_option
def ep(directory):
    """Run on the Epilepsy graph."""
    graph = get_ep_graph()
    run_epicom(graph, directory)


@main.command()
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
def multi(output):
    """Run on all graphs."""
    graphs = [
        get_ad_graph(),
        get_ep_graph(),
        get_pd_graph(),
    ]

    multi_run_epicom(graphs, output)


if __name__ == '__main__':
    main()
