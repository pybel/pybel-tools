# -*- coding: utf-8 -*-

"""

1. Download all physical protein-protein interactions from STRING
2. Filter a given network to all protein => protein interactions
3. Output all relations that don't appear

"""


def filter_ppi(graph):
    """Returns an iterator over directly increasing or directly decreasing relations between
    protein abundances in the given BEL Graph

    :param graph: A BEL Graph
    :type graph: pybel.BELGraph
    :return: A list of pairs of proteins that participate in direct protein-protein interactions
    :rtype: list
    """


def download_string(path):
    """Downloads the raw data from STRING

    :param path: The path to the file where to save the raw data
    :type path: str
    """
