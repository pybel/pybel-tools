# -*- coding: utf-8 -*-

"""
This module contains all of the services necessary through the PyBEL API Definition, backed by a network dictionary
"""

import logging
import time

from collections import Counter, defaultdict
from sqlalchemy import func

from pybel.manager.models import Network
from pybel.struct import union
from .mutation.metadata import add_canonical_names, enrich_pubmed_citations
from .utils import min_tanimoto_set_similarity

log = logging.getLogger(__name__)


class DatabaseService(object):
    """The dictionary service contains functions that implement the PyBEL API with a in-memory backend using
    dictionaries.
    """

    def __init__(self, manager, autocache=False):
        """
        :param pybel.manager.Manager manager: A cache manager
        """
        self.manager = manager

        #: dictionary of {int id: BELGraph graph}
        self.networks = {}

        #: A dictionary from {int network_id: {int network_id: float tanimoto node overlap}}
        self.overlap_cache = defaultdict(Counter)

        if autocache:
            self.cache_networks()

        log.info('initialized dictionary service')

    def clear(self):
        """Clear the cache and start over"""
        self.networks.clear()
        self.overlap_cache.clear()

    def _add_network(self, network_id, force_reload=False, eager=False, infer_origin=False):
        """Adds a network to the module-level cache from the underlying database

        :param int network_id: The identifier for the graph
        :param bool force_reload: Should the graphs be reloaded even if it has already been cached?
        :param bool eager: Should data be calculated/loaded eagerly?
        :param bool infer_origin: Should the central dogma be inferred for proteins, RNA, and miRNA?
        """
        if network_id in self.networks and not force_reload:
            return self.networks[network_id]

        log.debug('caching network [%s]', network_id)

        graph = self.manager.get_graph_by_id(network_id)

        t = time.time()

        log.info(
            'caching network [%s] %s v%s',
            network_id,
            graph.name,
            graph.version,
        )

        log.debug('adding canonical names')
        add_canonical_names(graph)

        if eager:
            log.debug('enriching citations')
            t = time.time()
            enrich_pubmed_citations(graph, manager=self.manager)
            log.debug('done enriching citations in %.2f seconds', time.time() - t)

        self.networks[network_id] = graph

        log.info(
            'cached (%d nodes, %d edges) in %.2f seconds',
            graph.number_of_nodes(),
            graph.number_of_edges(),
            time.time() - t
        )

        return graph

    def cache_networks(self, force_reload=False, eager=False, infer_origin=False):
        """This function needs to get all networks from the graph cache manager and make a dictionary

        :param bool force_reload: Should all graphs be reloaded even if they have already been cached?
        :param bool eager: Should difficult to preload features be calculated?
        :param bool infer_origin: Should the central dogma be inferred for proteins, RNA, and miRNA?
        """
        query = self.manager.session.query(Network).group_by(Network.name)

        if not force_reload and self.networks:
            query = query.filter(Network.id.notin_(self.networks))

        for network in query.having(func.max(Network.created)).order_by(Network.created.desc()).all():
            self._add_network(
                network.id,
                force_reload=force_reload,
                eager=eager,
                infer_origin=infer_origin,
            )

        log.info(
            'The cache has %d networks, %d nodes, %d edges',
            self.manager.count_networks(),
            self.manager.count_nodes(),
            self.manager.count_edges()
        )

    # Graph selection functions

    def get_graph_by_id(self, network_id):
        """Gets a network by its ID or super network if identifier is not specified

        :param int network_id: The internal ID of the network to get
        :rtype: pybel.BELGraph
        """
        return self._add_network(network_id)

    def get_graphs_by_ids(self, network_ids):
        """Gets a list of networks given the ids

        :param list[int] network_ids: A list of network identifiers
        :rtype: list[pybel.BELGraph]
        """
        return [
            self.get_graph_by_id(network_id)
            for network_id in network_ids
        ]

    def get_graph_by_ids(self, network_ids):
        """Gets a networks by a list of database identifiers

        :param list[int] network_ids: A list of network identifiers
        :rtype: pybel.BELGraph
        """
        if len(network_ids) == 1:
            return self.get_graph_by_id(network_ids[0])

        return union(self.get_graphs_by_ids(network_ids))

    def get_node_overlap(self, network_id):
        """Calculates overlaps to all other networks in the database

        :param int network_id: The network database identifier
        :return: A dictionary from {int network_id: float similarity} for this network to all other networks
        :rtype: collections.Counter[int,float]
        """
        t = time.time()

        network = self.manager.get_network_by_id(network_id)
        nodes = set(node.id for node in network.nodes)

        for other_network in self.manager.list_recent_networks():
            if other_network.id == network_id or other_network.id in self.overlap_cache[network_id]:
                continue

            other_network_nodes = set(node.id for node in other_network.nodes)

            overlap = min_tanimoto_set_similarity(nodes, other_network_nodes)

            self.overlap_cache[network_id][other_network.id] = overlap
            self.overlap_cache[other_network.id][network_id] = overlap

        log.debug('Cached node overlaps for network %s in %.2f seconds', network_id, time.time() - t)

        return self.overlap_cache[network_id]

    def forget_network(self, network_id):
        """Removes all cached data from the given network id"""
        if network_id in self.networks:
            del self.networks[network_id]

        if network_id in self.overlap_cache:
            del self.overlap_cache[network_id]

        for cached_network_id in self.overlap_cache:
            if network_id in self.overlap_cache[cached_network_id]:
                del self.overlap_cache[cached_network_id][network_id]
