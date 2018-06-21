# -*- coding: utf-8 -*-

"""

This module assists in running complex workflows on BEL graphs.

Example Pipeline #1
~~~~~~~~~~~~~~~~~~~
This example shows a pipeline that acquires a subgraph and finds possible additions to the subgraph.

>>> network = ...
>>> example = Pipeline()
>>> example.append('get_subgraph_by_annotation_value', 'Subgraph', 'Blood vessel dilation subgraph')
>>> example.append('enrich_unqualified')
>>> example.append('infer_central_dogma')
>>> example.append('expand_periphery')
>>> result = example.run(network)

Example Pipeline #2
~~~~~~~~~~~~~~~~~~~
This example shows how additional data can be integrated into a graph.

>>> from pybel.constants import PROTEIN
>>> network = ...
>>> example = Pipeline()
>>> example.append('infer_central_dogma')
>>> example.append('enrich_unqualified')
>>> example.append('infer_central_dogma')
>>> example.append('expand_periphery')
>>> example.append('expand_nodes_neighborhoods', [(PROTEIN, 'HGNC', 'AKT1'), (PROTEIN, 'HGNC', 'AKT2')])
>>> result = example.run(network)

Example Pipeline #3
~~~~~~~~~~~~~~~~~~~
This example shows how the results from multiple pipelines can be combine.

>>> network = ...
>>> pipeline_a = Pipeline()
>>> pipeline_a.append('get_subgraph_by_annotation_value', 'Subgraph', 'Blood vessel dilation subgraph')
>>> pipeline_b = Pipeline()
>>> pipeline_b.append('get_subgraph_by_annotation_value', 'Subgraph', 'Tau protein subgraph')
>>> pipeline_c = Pipeline.union(pipeline_a, pipeline_b)
>>> result = pipeline_c.run(network)

"""

import logging

from pybel.struct.pipeline.exc import MissingPipelineFunctionError
from pybel.struct.pipeline.pipeline import Pipeline, mapped

__all__ = [
    'splitter',
]

log = logging.getLogger(__name__)

#: A map of function names to functions that split graphs into dictionaries of graphs
splitter_map = {}


def assert_is_mapped_to_pipeline(name):
    """
    :param str name:
    :raises: MissingPipelineFunctionError
    """
    if name not in mapped:
        raise MissingPipelineFunctionError('{} is not registered as a pipeline function'.format(name))


def splitter(f):
    """A function decorator that signifies a function that takes in a graph and returns a dictionary of keys to graphs

    :param types.FunctionType f: A function
    :rtype: types.FunctionType
    """
    splitter_map[f.__name__] = f
    return f
