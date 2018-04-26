# -*- coding: utf-8 -*-

import os

import pybel

__all__ = [
    'get_ad_graph',
]


def get_bms_base():
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


def get_aetionomy_path():
    path = os.path.join(get_bms_base(), 'aetionomy')
    assert os.path.exists(path)
    return path


def get_ad_path():
    path = os.path.join(get_aetionomy_path(), 'alzheimers', 'alzheimers.bel')
    assert os.path.exists(path)
    return path


def get_ad_pickle_path():
    return os.path.join(get_aetionomy_path(), 'alzheimers', 'alzheimers.gpickle')


def get_ad_graph():
    pickle_path = get_ad_pickle_path()

    if os.path.exists(pickle_path):
        return pybel.from_pickle(pickle_path)

    raise RuntimeError('need to compile the AD knowledge assembly')


def get_pd_path():
    path = os.path.join(get_aetionomy_path(), 'parkinsons', 'parkinsons.bel')
    assert os.path.exists(path)
    return path


def get_pd_pickle_path():
    return os.path.join(get_aetionomy_path(), 'parkinsons', 'parkinsons.gpickle')


def get_pd_graph():
    pickle_path = get_pd_pickle_path()

    if os.path.exists(pickle_path):
        return pybel.from_pickle(pickle_path)

    raise RuntimeError('need to compile the PD knowledge assembly')


def get_ep_path():
    path = os.path.join(get_aetionomy_path(), 'epilepsy', 'epilepsy.bel')
    assert os.path.exists(path)
    return path


def get_ep_pickle_path():
    return os.path.join(get_aetionomy_path(), 'epilepsy', 'epilepsy.gpickle')


def get_ep_graph():
    pickle_path = get_ep_pickle_path()

    if os.path.exists(pickle_path):
        return pybel.from_pickle(pickle_path)

    raise RuntimeError('need to compile the Epilepsy knowledge assembly')
