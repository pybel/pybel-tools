# -*- coding: utf-8 -*-

"""This module contains functions that help overlay tabular data to nodes in a graph"""

import logging
from collections import defaultdict

import numpy as np

from pybel.constants import NAME
from pybel.struct.filters import filter_nodes
from .. import pipeline
from ..filters.node_filters import function_namespace_inclusion_builder

__all__ = [
    'overlay_data',
    'overlay_type_data',
]

log = logging.getLogger(__name__)


@pipeline.in_place_mutator
def overlay_data(graph, data, label, overwrite=False):
    """Overlays tabular data on the network

    :param pybel.BELGraph graph: A BEL Graph
    :param dict data: A dictionary of {tuple node: data for that node}
    :param str label: The annotation label to put in the node dictionary
    :param bool overwrite: Should old annotations be overwritten?
    """
    for node, value in data.items():
        if node not in graph:
            log.debug('%s not in graph', node)
            continue
        elif label in graph.node[node] and not overwrite:
            log.debug('%s already on %s', label, node)
            continue
        graph.node[node][label] = value


# TODO switch label to be kwarg with default value DATA_WEIGHT
@pipeline.in_place_mutator
def overlay_type_data(graph, data, label, func, namespace, overwrite=False, impute=None):
    """Overlays tabular data on the network for data that comes from an data set with identifiers that lack
    namespaces.

    For example, if you want to overlay differential gene expression data from a table, that table
    probably has HGNC identifiers, but no specific annotations that they are in the HGNC namespace or
    that the entities to which they refer are RNA.

    :param pybel.BELGraph graph: A BEL Graph
    :param dict[str,float] dict data: A dictionary of {name: data}
    :param str label: The annotation label to put in the node dictionary
    :param str func: The function of the keys in the data dictionary
    :param str namespace: The namespace of the keys in the data dictionary
    :param bool overwrite: Should old annotations be overwritten?
    :param Optional[float] impute: The value to use for missing data
    """
    new_data = {
        node: data.get(graph.node[node][NAME], impute)
        for node in filter_nodes(graph, function_namespace_inclusion_builder(func, namespace))
    }

    overlay_data(graph, new_data, label, overwrite=overwrite)


def load_differential_gene_expression(data_path, gene_symbol_column='Gene.symbol', logfc_column='logFC',
                                      aggregator=None):
    """Quick and dirty loader for differential gene expression data

    :param str data_path: The path to the CSV
    :param str gene_symbol_column: The header of the gene symbol column in the data frame
    :param str logfc_column: The header of the log-fold-change column in the data frame
    :param aggregator: A function that aggregates a list of differential gene expression values. Defaults to
                       :func:`numpy.median`. Could also use: :func:`numpy.mean`, :func:`numpy.average`,
                       :func:`numpy.min`, or :func:`numpy.max`
    :type aggregator: Optional[list[float] -> float]
    :return: A dictionary of {gene symbol: log fold change}
    :rtype: dict
    """
    import pandas as pd
    aggregator = aggregator or np.median

    df = pd.read_csv(data_path)
    df = df.loc[df[gene_symbol_column].notnull(), [gene_symbol_column, logfc_column]]

    values = defaultdict(list)

    for _, k, v in df.itertuples():
        values[k].append(v)

    return {
        k: aggregator(vs)
        for k, vs in values.items()
    }
