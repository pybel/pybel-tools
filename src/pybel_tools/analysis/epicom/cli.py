# -*- coding: utf-8 -*-

import os
import sys

import click

from .algorithm import multi_run_epicom, run_epicom
from ..neurommsig import get_ad_graph, get_ep_graph, get_pd_graph


@click.group()
def main():
    """Run EpiCom-Reloaded"""


@main.command()
@click.option('-d', '--directory', default=os.getcwd())
def ad(directory):
    """Run EpiCom on AD"""
    graph = get_ad_graph()
    run_epicom(graph, directory)


@main.command()
@click.option('-d', '--directory', default=os.getcwd())
def pd(directory):
    """Run EpiCom on PD"""
    graph = get_pd_graph()
    run_epicom(graph, directory)


@main.command()
@click.option('-d', '--directory', default=os.getcwd())
def ep(directory):
    """Run EpiCom on Epilepsy"""
    graph = get_ep_graph()
    run_epicom(graph, directory)


@main.command()
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
def multi(output):
    """Run EpiCom on all of NeuroMMSig"""
    graphs = [
        get_ad_graph(),
        get_ep_graph(),
        get_pd_graph(),
    ]

    multi_run_epicom(graphs, output)


if __name__ == '__main__':
    main()
