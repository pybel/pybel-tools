# -*- coding: utf-8 -*-

"""This module contains functions that help overlay tabular data to nodes in a graph"""

from collections import defaultdict

import logging
import numpy as np

from pybel.constants import NAME
from pybel.struct.filters import filter_nodes
from .. import pipeline
from ..constants import WEIGHT
from ..filters.node_filters import function_namespace_inclusion_builder

__all__ = [
    'overlay_data',
    'overlay_type_data',
    'load_differential_gene_expression',
]

log = logging.getLogger(__name__)


@pipeline.in_place_mutator
def overlay_data(graph, data, label=None, overwrite=False):
    """Overlays tabular data on the network

    :param pybel.BELGraph graph: A BEL Graph
    :param dict data: A dictionary of {tuple node: data for that node}
    :param Optional[str] label: The annotation label to put in the node dictionary
    :param bool overwrite: Should old annotations be overwritten?
    """
    if label is None:
        label = WEIGHT

    for node, value in data.items():
        if node not in graph:
            log.debug('%s not in graph', node)
            continue

        if label in graph.node[node] and not overwrite:
            log.debug('%s already on %s', label, node)
            continue

        graph.node[node][label] = value


@pipeline.in_place_mutator
def overlay_type_data(graph, data, func, namespace, label=None, overwrite=False, impute=None):
    """Overlay tabular data on the network for data that comes from an data set with identifiers that lack
    namespaces.

    For example, if you want to overlay differential gene expression data from a table, that table
    probably has HGNC identifiers, but no specific annotations that they are in the HGNC namespace or
    that the entities to which they refer are RNA.

    :param pybel.BELGraph graph: A BEL Graph
    :param dict[str,float] dict data: A dictionary of {name: data}
    :param str func: The function of the keys in the data dictionary
    :param str namespace: The namespace of the keys in the data dictionary
    :param Optional[str] label: The annotation label to put in the node dictionary
    :param bool overwrite: Should old annotations be overwritten?
    :param Optional[float] impute: The value to use for missing data
    """
    new_data = {
        node: data.get(graph.node[node][NAME], impute)
        for node in filter_nodes(graph, function_namespace_inclusion_builder(func, namespace))
    }

    overlay_data(graph, new_data, label=label, overwrite=overwrite)


def load_differential_gene_expression(path, gene_symbol_column='Gene.symbol', logfc_column='logFC', aggregator=None):
    """Load and preprocess a differential gene expression data.

    :param str path: The path to the CSV
    :param str gene_symbol_column: The header of the gene symbol column in the data frame
    :param str logfc_column: The header of the log-fold-change column in the data frame
    :param aggregator: A function that aggregates a list of differential gene expression values. Defaults to
                       :func:`numpy.median`. Could also use: :func:`numpy.mean`, :func:`numpy.average`,
                       :func:`numpy.min`, or :func:`numpy.max`
    :type aggregator: Optional[list[float] -> float]
    :return: A dictionary of {gene symbol: log fold change}
    :rtype: dict[str,float]
    """
    import pandas as pd

    if aggregator is None:
        aggregator = np.median

    # Load the data frame
    df = pd.read_csv(path)

    # Check the columns exist in the data frame
    assert gene_symbol_column in df.columns
    assert logfc_column in df.columns

    # throw away columns that don't have gene symbols - these represent control sequences
    df = df.loc[df[gene_symbol_column].notnull(), [gene_symbol_column, logfc_column]]

    values = defaultdict(list)

    for _, gene_symbol, log_fold_change in df.itertuples():
        values[gene_symbol].append(log_fold_change)

    return {
        gene_symbol: aggregator(log_fold_changes)
        for gene_symbol, log_fold_changes in values.items()
    }
