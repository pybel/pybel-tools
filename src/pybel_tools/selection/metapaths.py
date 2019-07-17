# -*- coding: utf-8 -*-

"""Metapath search utilities.

A metapath can be defined with two levels of granularity:

- Low: A list of BEL functions representing the types of entities in a given path
- High: An alternating list of BEL functions and BEL relations representing the types of entities in a given path and
  their relations

"""

from functools import lru_cache
from typing import Iterable, List, Tuple

from pybel import BELGraph
from pybel.dsl import BaseEntity

__all__ = [
    'convert_path_to_metapath',
    'yield_walks_exhaustive',
    'iterate_simple_metapaths',
]


def convert_path_to_metapath(nodes: Iterable[BaseEntity]) -> List[str]:
    """Convert a list of nodes to their corresponding functions.

    :param nodes: A list of BEL node tuples
    """
    return [
        node.function
        for node in nodes
    ]


@lru_cache(maxsize=None)
def yield_walks_exhaustive(graph, node, length: int) -> Iterable[Tuple]:
    """Get all walks under a given length starting at a given node.

    :param networkx.Graph graph: A graph
    :param node: Starting node
    :param length: The length of walks to get
    :return: A list of paths
    :rtype: list[tuple]
    """
    if 0 == length:
        return (node,),

    return (
        (node, relation) + path
        for neighbor in graph[node]
        for path in yield_walks_exhaustive(graph, neighbor, length - 1)
        if node not in path
        for relation in set(graph[node][neighbor][key]['relation'] for key in graph[node][neighbor])
    )


def iterate_simple_metapaths(
        graph: BELGraph,
        node: BaseEntity,
        simple_metapath: List[str],
) -> Iterable[Tuple[BaseEntity, ...]]:
    """Match a simple metapath starting at the given node.

    :param graph: A BEL graph
    :param node: A BEL node
    :param simple_metapath: A list of BEL Functions
    :return: An iterable over paths from the node matching the metapath

    .. code-block:: python

        from hbp_knowledge import get_graph
        from pybel.dsl import Protein
        from pybel.constants import PROTEIN

        graph = get_graph()

        mapt = Protein('HGNC', 'MAPT')

        for simple_metapath in iterate_simple_metapaths(graph, mapt, [PROTEIN, PROTEIN]):
            print(*(str(node) for node in simple_metapath))
    """
    if 0 == len(simple_metapath):
        yield node,

    else:
        for neighbor in graph[node]:
            if neighbor.function != simple_metapath[0]:
                continue
            for path in iterate_simple_metapaths(graph, neighbor, simple_metapath[1:]):
                if node not in path:
                    yield (node,) + path


def convert_simple_walk(simple_walk: Iterable[BaseEntity]) -> List[str]:
    """Convert a walk into a sequence of BEL functions.

    :param simple_walk: An iterable of BEL nodes
    :return: A list of BEL functions of the walk
    """
    return [
        node.function
        for node in simple_walk
    ]


def match_complex_metapath(graph, node, complex_metapath):
    """Match a complex metapath starting at the given node.

    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :param list[str] complex_metapath: An iterable of alternating BEL nodes and relations
    :return: An iterable over paths from the node matching the metapath
    :rtype: iter[tuple]
    """
    raise NotImplementedError


def convert_complex_walk(graph: BELGraph, complex_walk):
    """Convert a walk into an alternative sequence of BEL functions and BEL relations.

     This result is starting and ending with a BEL function

    :param graph: A BEL graph
    :param iter[tuple] complex_walk: An iterable of alternating BEL nodes and relations
    :return: An alternating list of BEL functions and relations of the walk
    :rtype: list[str]
    """
    raise NotImplementedError


def _main(max_length=3):
    from hbp_knowledge import get_graph

    graph = get_graph()
    graph.summarize()

    c = 0
    for node in graph:
        walks = yield_walks_exhaustive(graph, node, max_length)
        for i, walk in enumerate(walks):
            print(i, *[str(n) for n in walk])
            c += 1

    print(f'Found {c} walks of max_length={max_length} on {graph}')


if __name__ == '__main__':
    _main()
