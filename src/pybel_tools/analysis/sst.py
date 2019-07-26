# -*- coding: utf-8 -*-

"""An implementation of the sampling of spanning trees (SST) algorithm [﻿Vasilyev2014]_.

.. [﻿Vasilyev2014] ﻿Vasilyev, D. M., *et al* (2014). `An algorithm for score aggregation over causal biological networks
                  based on random walk sampling. <https://doi.org/10.1186/1756-0500-7-516>`_ BMC Research Notes, 7, 516.
"""


def get_random_walk_spanning_tree(graph):
    """Generate a spanning tree from the directed graph using the random walk approach.

    This was proposed independently by by Broder (1989) and Aldous (1990). It simply generates random walks until all
    nodes have been covered.

    Algorithm:

    1. Choose a starting vertex s arbitrarily. Set T_V ← {s} and T_E ← ∅.
    2. Do a simple random walk starting at s. Whenever we cross an edge e = {u, v} with v ∈ V ,
        add v to TV and add e to TE.
    3. Stop the random walk when TV = V . Output T = (T_V , T_E) as our spanning tree

    :param networkx.DiGraph graph: The input graph
    :rtype: networkx.DiGraph

    .. seealso::

        - https://math.dartmouth.edu/~pw/math100w13/kothari.pdf
        - http://keyulux.com/pdf/spanning_tree.pdf

    """
    raise NotImplementedError
