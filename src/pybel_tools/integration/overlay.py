# -*- coding: utf-8 -*-

"""This module contains functions that help overlay tabular data to nodes in a graph"""

import logging
from collections import defaultdict
from typing import Any, Callable, List, Mapping, Optional

import numpy as np
import pandas as pd

from pybel import BELGraph
from pybel.constants import NAME
from pybel.dsl import BaseEntity
from pybel.struct.filters import filter_nodes
from pybel.struct.pipeline import in_place_transformation
from ..constants import WEIGHT
from ..filters.node_filters import function_namespace_inclusion_builder

__all__ = [
    'overlay_data',
    'overlay_type_data',
    'load_differential_gene_expression',
]

log = logging.getLogger(__name__)


@in_place_transformation
def overlay_data(graph: BELGraph,
                 data: Mapping[BaseEntity, Any],
                 label: Optional[str] = None,
                 overwrite: bool = False,
                 ) -> None:
    """Overlays tabular data on the network

    :param graph: A BEL Graph
    :param data: A dictionary of {tuple node: data for that node}
    :param label: The annotation label to put in the node dictionary
    :param overwrite: Should old annotations be overwritten?
    """
    if label is None:
        label = WEIGHT

    for node, value in data.items():
        if node not in graph:
            log.debug('%s not in graph', node)
            continue

        if label in graph.nodes[node] and not overwrite:
            log.debug('%s already on %s', label, node)
            continue

        graph.nodes[node][label] = value


@in_place_transformation
def overlay_type_data(graph: BELGraph,
                      data: Mapping[str, float],
                      func: str,
                      namespace: str,
                      label: Optional[str] = None,
                      overwrite: bool = False,
                      impute: Optional[float] = None,
                      ) -> None:
    """Overlay tabular data on the network for data that comes from an data set with identifiers that lack
    namespaces.

    For example, if you want to overlay differential gene expression data from a table, that table
    probably has HGNC identifiers, but no specific annotations that they are in the HGNC namespace or
    that the entities to which they refer are RNA.

    :param graph: A BEL Graph
    :param dict data: A dictionary of {name: data}
    :param func: The function of the keys in the data dictionary
    :param namespace: The namespace of the keys in the data dictionary
    :param label: The annotation label to put in the node dictionary
    :param overwrite: Should old annotations be overwritten?
    :param impute: The value to use for missing data
    """
    new_data = {
        node: data.get(node[NAME], impute)
        for node in filter_nodes(graph, function_namespace_inclusion_builder(func, namespace))
    }

    overlay_data(graph, new_data, label=label, overwrite=overwrite)


def load_differential_gene_expression(path: str,
                                      gene_symbol_column: str = 'Gene.symbol',
                                      logfc_column: str = 'logFC',
                                      aggregator: Optional[Callable[[List[float]], float]] = None,
                                      ) -> Mapping[str, float]:
    """Load and pre-process a differential gene expression data.

    :param path: The path to the CSV
    :param gene_symbol_column: The header of the gene symbol column in the data frame
    :param logfc_column: The header of the log-fold-change column in the data frame
    :param aggregator: A function that aggregates a list of differential gene expression values. Defaults to
                       :func:`numpy.median`. Could also use: :func:`numpy.mean`, :func:`numpy.average`,
                       :func:`numpy.min`, or :func:`numpy.max`
    :return: A dictionary of {gene symbol: log fold change}
    """
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
