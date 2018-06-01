# -*- coding: utf-8 -*-

"""This module contains the CLI to export NeuroMMSig as BEL.

To run, type :code:`python3 -m pybel_tools.analysis.neurommsig.cli` in the command line
"""

import logging
import os

import click

from .export import get_nift_values, mesh_alzheimer, mesh_parkinson, preprocess, write_neurommsig_bel
from .utils import get_bms_base, get_neurommsig_base

log = logging.getLogger(__name__)


@click.command()
def main():
    logging.basicConfig(level=logging.INFO)
    log.setLevel(logging.INFO)

    bms_base = get_bms_base()
    neurommsig_base = get_neurommsig_base()
    neurommsig_excel_dir = os.path.join(neurommsig_base, 'resources', 'excels', 'neurommsig')

    nift_values = get_nift_values()

    log.info('Starting Alzheimers')

    ad_path = os.path.join(neurommsig_excel_dir, 'alzheimers', 'alzheimers.xlsx')
    ad_df = preprocess(ad_path)
    with open(os.path.join(bms_base, 'aetionomy', 'alzheimers', 'neurommsigdb_ad.bel'), 'w') as ad_file:
        write_neurommsig_bel(ad_file, ad_df, mesh_alzheimer, nift_values)

    log.info('Starting Parkinsons')

    pd_path = os.path.join(neurommsig_excel_dir, 'parkinsons', 'parkinsons.xlsx')
    pd_df = preprocess(pd_path)
    with open(os.path.join(bms_base, 'aetionomy', 'parkinsons', 'neurommsigdb_pd.bel'), 'w') as pd_file:
        write_neurommsig_bel(pd_file, pd_df, mesh_parkinson, nift_values)


if __name__ == '__main__':
    main()
