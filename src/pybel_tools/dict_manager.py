# -*- coding: utf-8 -*-

"""This class contains an alternate implementation of the PyBEL database manager that only stores graphs in memory."""

from typing import Iterable, List, Optional

from pybel import BELGraph, Manager
from pybel.manager.models import Network


class _Namespace:
    pass


class DictManager(Manager):
    """A dictionary-based implementation of the PyBEL Manager."""

    def __init__(self, connection: Optional[str] = None):
        super(DictManager, self).__init__(connection=connection)

        self.universe = None
        self.networks = {}
        self.disease_to_id = {}
        self.hash_to_node = {}

    def insert_graph(self, graph: BELGraph, **_kwargs) -> Network:
        """Insert a graph and return the resulting ORM object (mocked)."""
        result = _Namespace()
        result.id = len(self.networks)

        self.networks[result.id] = graph

        return result

    def get_graph_by_id(self, network_id: int) -> BELGraph:
        """Get a graph by its identifier."""
        return self.networks[network_id]

    def get_graphs_by_ids(self, network_ids: Iterable[int]) -> List[BELGraph]:
        """Get several graphs by their identifiers."""
        return [
            self.networks[network_id]
            for network_id in network_ids
        ]
