# -*- coding: utf-8 -*-

import os
import sys

import click

import pybel
from pybel_tools.analysis.epicom import run_epicom

bms_base = os.environ['BMS_BASE']


@click.command()
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
def main(output):
    ad_pickle_path = os.path.join(bms_base, 'aetionomy', 'alzheimers', 'alzheimers.gpickle')
    assert os.path.exists(ad_pickle_path)

    ep_pickle_path = os.path.join(bms_base, 'aetionomy', 'epilepsy', 'epilepsy.gpickle')
    assert os.path.exists(ep_pickle_path)

    graphs = [
        pybel.from_pickle(ad_pickle_path),
        pybel.from_pickle(ep_pickle_path),
    ]

    run_epicom(graphs, output)


if __name__ == '__main__':
    main()
