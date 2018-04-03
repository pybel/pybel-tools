# -*- coding: utf-8 -*-

from collections import defaultdict

import networkx as nx
import pandas as pd

__all__ = [
    'pseudotoposort_runner',
]


# reimplement the topological sort but be permissive of non-DAGS
def _topological_sort_recursive(G, reverse=False):
    """Return a list of nodes in topological sort order.

    A topological sort is a nonunique permutation of the nodes such
    that an edge from u to v implies that u appears before v in the
    topological sort order.

    Parameters
    ----------
    G : NetworkX digraph

    nbunch : container of nodes (optional)
       Explore graph in specified order given in nbunch

    reverse : bool, optional
        Return postorder instead of preorder if True.
        Reverse mode is a bit more efficient.

    Raises
    ------
    NetworkXError
       Topological sort is defined for directed graphs only. If the
       graph G is undirected, a NetworkXError is raised.

    NetworkXUnfeasible
        If G is not a directed acyclic graph (DAG) no topological sort
        exists and a NetworkXUnfeasible exception is raised.

    Notes
    -----
    This is a recursive version of topological sort.

    See also
    --------
    topological_sort
    is_directed_acyclic_graph

    """
    if not G.is_directed():
        raise nx.NetworkXError(
            "Topological sort not defined on undirected graphs.")

    def _dfs(n):
        ancestors.add(n)

        for w in G[n]:
            if w not in explored:
                _dfs(w)

        ancestors.remove(n)
        explored.add(n)
        order.append(n)

    ancestors = set()
    explored = set()
    order = []

    for v in list(G):
        if v not in explored:
            _dfs(v)

    if reverse:
        return order
    else:
        return list(reversed(order))


def pseudotoposort_runner(graph, trials=None):
    """

    1. Randomly select a spanning DAG with a randomized DFS
    2. Return topological sort
    3. Repeat multiple times

    Return a matrix with one axis - nodes, one axis - position

    :param graph:
    :param Optional[int] trials: Number of trials to run. Defaults to 100.
    :rtype: list[list]
    """
    if nx.is_directed_acyclic_graph(graph):
        raise ValueError('Graph is a dag!')

    r = defaultdict(int)

    for _ in range(trials or 100):
        for i, node in enumerate(nx.dfs_postorder_nodes(graph)):
            r[i, node] += 1

    return pd.DataFrame(r)
