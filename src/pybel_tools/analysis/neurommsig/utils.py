# -*- coding: utf-8 -*-

"""Utilities for getting NeuroMMSig graphs."""

import os

from bel_repository import BELRepository
from pybel import BELGraph

__all__ = [
    'get_ad_graph',
    'get_pd_graph',
    'get_ep_graph',
    'get_bms_base',
    'get_neurommsig_base',
]


def get_ad_graph(*args, **kwargs) -> BELGraph:
    repo = BELRepository(os.path.join(get_aetionomy_directory(), 'alzheimers'))
    return repo.get_graph(*args, **kwargs)


def get_pd_graph(*args, **kwargs) -> BELGraph:
    repo = BELRepository(os.path.join(get_aetionomy_directory(), 'parkinsons'))
    return repo.get_graph(*args, **kwargs)


def get_ep_graph(*args, **kwargs) -> BELGraph:
    repo = BELRepository(os.path.join(get_aetionomy_directory(), 'epilepsy'))
    return repo.get_graph(*args, **kwargs)


def get_aetionomy_directory() -> str:
    aetionomy_bel_directory = os.environ.get('AETIONOMY_BEL_DIRECTORY')
    if not aetionomy_bel_directory:
        aetionomy_bel_directory = get_bms_base()

    path = os.path.join(aetionomy_bel_directory, 'aetionomy')
    assert os.path.exists(path)
    return path


def get_bms_base() -> str:
    """Get the path to the BMS git repository from the environment or thrown an exception.

    :raises: RuntimeError
    """
    bms_base = os.environ.get('BMS_BASE')

    if bms_base is None:
        raise RuntimeError("""
        
        BMS_BASE environment variable is not set.
        
        1. Change to home directory with `cd ~`
        2. Create a folder called dev with `mkdir dev` and enter with `cd dev`
        3. Clone BMS with: `git clone https://gitlab.scai.fraunhofer.de/bio/bms.git`
        4. Set the environment variable BMS_BASE in the current session with `export BMS_BASE="~/dev/bms"` or add
           to your .bashrc
        
        """)

    return bms_base


def get_neurommsig_base() -> str:
    """Get the path to the NeuroMMSig git repository from the environment or thrown an exception.

    :raises: RuntimeError
    """
    neurommsig_base = os.environ.get('NEUROMMSIG_BASE')

    if neurommsig_base is None:
        raise RuntimeError("""
        
        NEUROMMSIG_BASE environment variable is not set.
        
        1. Change to home directory with `cd ~`
        2. Create a folder called dev with `mkdir dev` and enter with `cd dev`
        3. Clone NeuroMMSig with `git clone https://gitlab.scai.fraunhofer.de/daniel.domingo.fernandez/neurommsig.git`
        4. Set the environment variable NEUROMMSIG_BASE in the current session with 
           `export NEUROMMSIG_BASE="~/dev/bms"` or add to your .bashrc
           
        """)

    return neurommsig_base
