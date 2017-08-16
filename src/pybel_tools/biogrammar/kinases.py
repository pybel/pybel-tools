# -*- coding: utf-8 -*-

"""

1. Download all enzyme classification data from BRENDA
2. Search for all edges with either the subject or object having a molecular activity
3. Check the corresponding nodes to have enzyme classifications matching the molecular activities

For example, kinase activities would be important

"""


def download_brenda(path):
    """Downloads the BRENDA data

    :param path: The path to the file where to save the raw data
    :type path: str
    """


def extract_molecular_activities(graph):
    """Extracts nodes that have molecular activities, such as kinase activities


    :param graph: A BEL graph
    :type graph: pybel.BELGraph
    :return: A list of nodes that have molecular activities
    :rtype: list
    """