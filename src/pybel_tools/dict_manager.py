# -*- coding: utf-8 -*-

"""This class contains an alternate implementation of the PyBEL database manager that only stores graphs
in memory
"""

from pybel.manager import Manager


class _Namespace(object):
    pass


class DictManager(Manager):
    """A dictionary-based implementation of the PyBEL Manager"""

    def __init__(self, connection=None):
        """
        :param Optional[str] connection:
        """
        super(DictManager, self).__init__(connection=connection)

        self.universe = None
        self.networks = {}
        self.disease_to_id = {}
        self.hash_to_node = {}

    def insert_graph(self, graph, **kwargs):
        """
        :param pybel.BELGraph graph:
        :param kwargs: Swallowed. Just used to match signature of pybel.Manager
        :rtype: Network
        """
        result = _Namespace()
        result.id = len(self.networks)

        self.networks[result.id] = graph

        return result

    def get_graph_by_id(self, network_id):
        """Returns a graph by its id

        :param int network_id:
        :rtype: pybel.BELGraph
        """
        return self.networks[network_id]

    def get_graphs_by_ids(self, network_ids):
        """
        :param list[int] network_ids:
        :rtype: list[pybel.BELGraph]
        """
        return [
            self.networks[network_id]
            for network_id in network_ids
        ]
