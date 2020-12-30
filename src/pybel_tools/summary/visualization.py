# -*- coding: utf-8 -*-

"""Functions for summarizing graphs.

This module contains functions that provide aggregate summaries of graphs including visualization with matplotlib,
printing summary information, and exporting summarized graphs
"""

import logging

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from pybel import BELGraph

__all__ = [
    'plot_summary_axes',
    'plot_summary',
]

logger = logging.getLogger(__name__)


def plot_summary_axes(graph: BELGraph, lax, rax, logx: bool = True):
    """Plot the graph summary statistics on the given axes.

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
    function_counter = graph.count.functions()
    relation_counter = graph.count.relations()

    function_df = pd.DataFrame.from_dict(dict(function_counter), orient='index').reset_index()
    function_df.columns = ['Function', 'Count']
    function_df.sort_values('Count', ascending=False, inplace=True)
    relation_df = pd.DataFrame.from_dict(dict(relation_counter), orient='index').reset_index()
    relation_df.columns = ['Relation', 'Count']
    relation_df.sort_values('Count', ascending=False, inplace=True)

    sns.barplot(x='Count', y='Function', data=function_df, ax=lax, orient='h')
    lax.set_title('Number of nodes: {}'.format(graph.number_of_nodes()))

    sns.barplot(x='Count', y='Relation', data=relation_df, ax=rax, orient='h')
    rax.set_title('Number of edges: {}'.format(graph.number_of_edges()))

    if logx:
        lax.set_xscale('log')
        rax.set_xscale('log')


def plot_summary(graph: BELGraph, logx: bool = True, **kwargs):
    """Plot your graph summary statistics.

    This function is a thin wrapper around :func:`plot_summary_axis`. It
    automatically takes care of building figures given matplotlib's pyplot module as an argument. After, you need
    to run :func:`plt.show`.

    :code:`plt` is given as an argument to avoid needing matplotlib as a dependency for this function

    Shows:

    1. Count of nodes, grouped by function type
    2. Count of edges, grouped by relation type

    :param kwargs: keyword arguments to give to :func:`plt.subplots`

    Example usage:

    >>> import matplotlib.pyplot as plt
    >>> from pybel import from_pickle
    >>> from pybel_tools.summary import plot_summary
    >>> graph = from_pickle('~/dev/bms/aetionomy/parkinsons.gpickle')
    >>> plot_summary(graph, figsize=(10, 4))
    >>> plt.show()
    """
    fig, (lax, rax) = plt.subplots(1, 2, **kwargs)

    plot_summary_axes(graph, lax, rax, logx=logx)
    plt.tight_layout()

    return fig, (lax, rax)
