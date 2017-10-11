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
from .mutation.metadata import add_canonical_names, add_identifiers, enrich_pubmed_citations
from .utils import (
    min_tanimoto_set_similarity,
)

log = logging.getLogger(__name__)


class DatabaseServiceBase:
    def __init__(self, manager):
        """
        :param pybel.manager.Manager manager: A cache manager
        """
        self.manager = manager


# TODO delete?
class QueryService(DatabaseServiceBase):
    def get_edges_by_network_id(self, network_id, offset_start=0, offset_end=500):
        network = self.manager.get_network_by_id(network_id)

        if offset_start > 0:
            offset_start -= 1

        if offset_start == offset_end:
            offset_end += 1

        edges = [
            edge.to_json(include_id=True)
            for edge in network.edges[offset_start:offset_end]
        ]

        result = {
            'network': network.to_json(include_id=True),
            'edges': edges,
            'number_of_edges': len(edges),
            'offset': {
                'start': offset_start,
                'end': offset_end
            }
        }

        return result

    def query_edges(self, network_id=None, pmid=None, statement=None, source=None, target=None, relation=None,
                    author=None, citation=None, annotations=None, offset_start=None, offset_end=None):
        """Provides a list of edges (nanopubs) filtered by the given parameters.

        :param int network_id: Primary identifier of the network in the PyBEL database. This can be obtained with the
                           get_networks function.
        :param str pmid: Pubmed identifier of a specific citation.
        :param str statement: A BEL statement that represents the needed edge.
        :param str source: A BEL term that describes the SUBJECT of the seeked relation.
        :param str target: A BEL term that describes the OBJECT of the seeked relation.
        :param str relation: The relation that is used in the seeked relationship.
        :param str author: An author that participated to the cited article.
        :param str or pybel.models.Citation citation: A citation that is used to back up the given relationship.
        :param dict annotations: A dictionary that describes an annotation that is the context of the relationship.
        :param int offset_start: The starting point of the offset (position in database). Defaults to 0
        :param int offset_end: The end point of the offset (position in database). Defaults to 500
        :rtype: list[]
        """
        offset_start = offset_start if offset_start is not None else 0
        offset_end = offset_end if offset_end is not None else 500

        if network_id:
            return self.get_edges_by_network_id(network_id, offset_start, offset_end)

        if author and citation is None and pmid is None:
            citation = self.manager.query_citations(author=author)
        elif pmid:
            citation = str(pmid)

        edges = self.manager.query_edges(
            bel=statement,
            source=source,
            target=target,
            relation=relation,
            citation=citation,
            annotation=annotations
        )

        result = {
            'number_of_edges': len(edges),
            'edges': [
                edge.to_json(include_id=True)
                for edge in edges
            ],
        }

        return result

    def query_nodes(self, network_id=None, node_id=None, bel=None, namespace=None, name=None, offset_start=0,
                    offset_end=500):
        """Provides a list of nodes filtered by the given parameters.

        :param network_id: Primary identifier of the network in the PyBEL database. This can be obtained with the
                           :func:`get_networks` function.
        :type network_id: int
        :param node_id: A node's database identifier
        :type node_id: int
        :param bel: A BEL term that describes the seeked node.
        :type bel: str
        :param namespace: A namespace keyword (e.g. HGNC) that represents a namespace that describes the nodes name.
        :type namespace: str
        :param name: The name of the node (biological entity).
        :type name: str
        :param offset_start: The starting point of the offset (position in database)
        :type offset_start: int
        :param offset_end: The end point of the offset (position in database)
        :type offset_end: int
        :return:
        :rtype:
        """
        if node_id:
            result = self.manager.query_nodes(node_id=node_id)
        elif network_id:
            network = self.manager.get_network_by_id(network_id)
            result = [node.data for node in network.nodes[offset_start:offset_end]]
        else:
            result = self.manager.query_nodes(bel=bel, namespace=namespace, name=name, as_dict_list=True)

        return result


class DatabaseService(QueryService):
    """The dictionary service contains functions that implement the PyBEL API with a in-memory backend using
    dictionaries.
    """

    def __init__(self, manager, autocache=False):
        """
        :param pybel.manager.Manager manager: A cache manager
        """
        super(DatabaseService, self).__init__(manager)

        #: dictionary of {int id: BELGraph graph}
        self.networks = {}

        #: A dictionary from {int network_id: {int network_id: float tanimoto node overlap}}
        self.overlap_cache = defaultdict(dict)

        if autocache:
            self.cache_networks()

        log.info('initialized dictionary service')

    def clear(self):
        """Clear the cache and start over"""
        self.__init__(self.manager)

    def _update_indexes(self, graph):
        """Updates identifiers for nodes based on addition order

        :param pybel.BELGraph graph: A BEL Graph
        """
        add_identifiers(graph)  # adds stable identifiers to nodes and edges

        for node, data in graph.nodes_iter(data=True):
            if isinstance(node, int):
                log.warning('nodes already converted to ids')
                return

    def _add_network(self, network_id, force_reload=False, eager=False, infer_origin=False):
        """Adds a network to the module-level cache from the underlying database

        :param int network_id: The identifier for the graph
        :param bool force_reload: Should the graphs be reloaded even if it has already been cached?
        :param bool eager: Should data be calculated/loaded eagerly?
        :param bool infer_origin: Should the central dogma be inferred for proteins, RNA, and miRNA?
        """
        log.debug('caching network [%s]', network_id)

        if network_id in self.networks and not force_reload:
            log.info('tried re-adding graph [%s]', network_id)
            return self.networks[network_id]

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

        log.debug('updating graph node/edge indexes')
        self._update_indexes(graph)

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
            'The database has %d networks, %s nodes, %s edges',
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
        if network_id not in self.networks:
            self._add_network(network_id)

        return self.networks[network_id]

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

        other_network_ids = {
            other_network.id
            for other_network in self.manager.list_recent_networks()
            if other_network.id not in self.overlap_cache[network_id]
        }

        if not other_network_ids:
            return Counter(self.overlap_cache[network_id])

        graph = self.get_graph_by_id(network_id)
        nodes = set(graph)

        for other_network_id in other_network_ids:
            if other_network_id == network_id:
                continue

            if other_network_id in self.overlap_cache[network_id]:
                continue

            other_network = self.get_graph_by_id(other_network_id)
            other_network_nodes = set(other_network)

            overlap = min_tanimoto_set_similarity(nodes, other_network_nodes)

            self.overlap_cache[network_id][other_network_id] = overlap
            self.overlap_cache[other_network_id][network_id] = overlap

        log.debug('Cached node overlaps for %s in %.2f seconds', graph, time.time() - t)

        return Counter(self.overlap_cache[network_id])

    def forget_network(self, network_id):
        """Removes all cached data from the given network id"""
        if network_id in self.networks:
            del self.networks[network_id]

        if network_id in self.overlap_cache:
            del self.overlap_cache[network_id]

        for cached_network_id in self.overlap_cache:
            if network_id in self.overlap_cache[cached_network_id]:
                del self.overlap_cache[cached_network_id][network_id]
