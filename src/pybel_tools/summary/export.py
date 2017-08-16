# -*- coding: utf-8 -*-

"""

This module contains functions that provide aggregate summaries of graphs including visualization with matplotlib,
printing summary information, and exporting summarized graphs

"""

from __future__ import print_function

import logging

import networkx as nx
import pandas as pd

from .edge_summary import count_relations
from .node_summary import count_functions
from .provenance import count_unique_citations, count_unique_authors

__all__ = [
    'plot_summary_axes',
    'plot_summary',
    'info_list',
    'info_str',
    'info_json',
    'print_summary',
]

log = logging.getLogger(__name__)


def plot_summary_axes(graph, lax, rax):
    """Plots your graph summary statistics on the given axes.

    After, you should run :func:`plt.tight_layout` and you must run :func:`plt.show` to view.

    Shows:
    1. Count of nodes, grouped by function type
    2. Count of edges, grouped by relation type

    :param pybel.BELGraph graph: A BEL graph
    :param lax: An axis object from matplotlib
    :param rax: An axis object from matplotlib

    Example usage:

    >>> import matplotlib.pyplot as plt
    >>> from pybel import from_pickle
    >>> from pybel_tools.summary import plot_summary_axes
    >>> graph = from_pickle('~/dev/bms/aetionomy/parkinsons.gpickle')
    >>> fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    >>> plot_summary_axes(graph, axes[0], axes[1])
    >>> plt.tight_layout()
    >>> plt.show()
    """
    ntc = count_functions(graph)
    etc = count_relations(graph)

    df = pd.DataFrame.from_dict(dict(ntc), orient='index')
    df_ec = pd.DataFrame.from_dict(dict(etc), orient='index')

    df.sort_values(0, ascending=True).plot(kind='barh', logx=True, ax=lax)
    lax.set_title('Number of nodes: {}'.format(graph.number_of_nodes()))

    df_ec.sort_values(0, ascending=True).plot(kind='barh', logx=True, ax=rax)
    rax.set_title('Number of edges: {}'.format(graph.number_of_edges()))


def plot_summary(graph, plt, **kwargs):
    """Plots your graph summary statistics. This function is a thin wrapper around :func:`plot_summary_axis`. It
    automatically takes care of building figures given matplotlib's pyplot module as an argument. After, you need
    to run :func:`plt.show`.

    :code:`plt` is given as an argument to avoid needing matplotlib as a dependency for this function

    Shows:

    1. Count of nodes, grouped by function type
    2. Count of edges, grouped by relation type

    :param pybel.BELGraph graph: A BEL graph
    :param plt: Give :code:`matplotlib.pyplot` to this parameter
    :param kwargs: keyword arguments to give to :func:`plt.subplots`

    Example usage:

    >>> import matplotlib.pyplot as plt
    >>> from pybel import from_pickle
    >>> from pybel_tools.summary import plot_summary
    >>> graph = from_pickle('~/dev/bms/aetionomy/parkinsons.gpickle')
    >>> plot_summary(graph, plt, figsize=(10, 4))
    >>> plt.show()
    """
    fig, axes = plt.subplots(1, 2, **kwargs)
    lax = axes[0]
    rax = axes[1]

    plot_summary_axes(graph, lax, rax)
    plt.tight_layout()

    return fig, axes


def info_list(graph):
    """Returns useful information about the graph as a list of tuples

    :param pybel.BELGraph graph: A BEL graph
    :rtype: list
    """
    number_nodes = graph.number_of_nodes()
    result = [
        ('Nodes', number_nodes),
        ('Edges', graph.number_of_edges()),
        ('Citations', count_unique_citations(graph)),
        ('Authors', count_unique_authors(graph)),
        ('Network density', nx.density(graph)),
        ('Components', nx.number_weakly_connected_components(graph)),
    ]

    try:
        result.append(('Average degree', sum(graph.in_degree().values()) / float(number_nodes)))
    except ZeroDivisionError:
        log.info('Graph has no nodes')

    if graph.warnings:
        result.append(('Compilation warnings', len(graph.warnings)))

    return result


def info_json(graph):
    """Returns useful information about the graph as a dictionary

    :param pybel.BELGraph graph: A BEL graph
    :rtype: dict
    """
    return dict(info_list(graph))


def info_str(graph):
    """Puts useful information about the graph in a string

    :param pybel.BELGraph graph: A BEL graph
    :rtype: str
    """
    return '\n'.join('{}: {}'.format(k, v) for k, v in info_list(graph))


def print_summary(graph, file=None):
    """Prints useful information about the graph

    :param pybel.BELGraph graph: A BEL graph
    :param file: A writeable file or file-like object. If None, defaults to :data:`sys.stdout`
    """
    print(info_str(graph), file=file)
