# -*- coding: utf-8 -*-

"""
This module contains all of the services necessary through the PyBEL API Definition, backed by a network dictionary
"""

import logging
import time

import networkx as nx
from collections import Counter, defaultdict
from functools import lru_cache
from sqlalchemy import func

from pybel import BELGraph
from pybel.canonicalize import node_to_bel, calculate_canonical_name
from pybel.constants import ID
from pybel.manager.models import Network, Annotation
from pybel.struct import left_full_join, union
from pybel.utils import edge_to_tuple
from .constants import CNAME
from .mutation.inference import infer_central_dogma
from .mutation.metadata import (
    parse_authors,
    add_canonical_names,
    enrich_pubmed_citations,
    add_identifiers,
)
from .summary.edge_summary import (
    count_pathologies,
    get_tree_annotations,
    get_annotations_containing_keyword
)
from .summary.provenance import (
    get_authors,
    get_pmid_by_keyword,
    get_authors_by_keyword,
    get_pubmed_identifiers
)
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


class QueryService(DatabaseServiceBase):
    def query_namespaces(self, network_id=None, offset_start=0, offset_end=500, name_list=False, keyword=None):
        """Provides a list of namespaces filtered by the given parameters.

        :param int network_id: Primary identifier of the network in the PyBEL database.
        :param int offset_start: The starting point of the offset (position in database)
        :param int offset_end: The end point of the offset (position in database)
        :param bool name_list: Flag that identifies if a list of all names in a namespace should be created.
        :param str keyword: The keyword used to identify a specific namespace. This is only used if name_list is True.
        :return: List of dictionaries that describe all namespaces.
        """
        result = []

        if network_id:
            network = self.manager.get_network_by_id(network_id)
            result = [
                namespace.to_json()
                for namespace in network.namespaces[offset_start:offset_end]
            ]

        if name_list and keyword:
            keyword_url_dict = {
                namespace.keyword: namespace.url
                for namespace in self.manager.list_namespaces()
            }
            namespace_url = keyword_url_dict[keyword]
            self.manager.ensure_namespace(url=namespace_url)

            definition = self.manager.get_namespace_by_url(namespace_url)
            namespace_data = definition.to_json()

            result = {
                'namespace_definition': namespace_data,
                'names': self.manager.namespace_cache[namespace_url]
            }

        if not result:
            result = [
                definition.to_json()
                for definition in self.manager.list_namespaces()
            ]

        return result

    def query_annotations(self, network_id=None, offset_start=0, offset_end=500, name_list=False, keyword=None):
        """Provides a list of annotations filtered by the given parameters.

        :param int network_id: Primary identifier of the network in the PyBEL database. This can be obtained with the
                           get_networks function.
        :param int offset_start: The starting point of the offset (position in database)
        :param int offset_end: The end point of the offset (position in database)
        :param bool name_list: Flag that identifies if a list of all names in a namespace should be created.
        :param str keyword: The keyword used to identify a specific namespace. This is only used if name_list is True.
        :return: List of dictionaries that describe all namespaces.
        """
        result = []

        if network_id:
            network = self.manager.get_network_by_id(network_id)
            result = [annotation.data for annotation in network.annotations[offset_start:offset_end]]

        if name_list:
            if name_list and keyword:
                keyword_url_dict = dict(self.manager.session.query(Annotation.keyword, Annotation.url).all())
                url = keyword_url_dict[keyword]
                self.manager.ensure_annotation(url=url)
                annotation = self.manager.session.query(Annotation).filter_by(url=url).one_or_none()
                annotation_data = annotation.to_json()

                result = {
                    'annotation_definition': annotation_data,
                    'annotations': self.manager.annotation_cache[url]
                }

        if len(result) == 0:
            result = [
                definition.to_json()
                for definition in self.manager.session.query(Annotation).all()
            ]

        return result

    def query_citations(self, network_id=None, author=None, offset_start=0, offset_end=500):
        """

        :param int network_id: Database identifier of the network in the PyBEL database
        :param str author: The name of an author that participated in creation of the citation.
        :param int offset_start: The starting point of the offset (position in database)
        :param int offset_end: The end point of the offset (position in database)
        :return:
        """
        result = []

        if network_id:
            network = self.manager.get_network_by_id(network_id)
            result = [
                citation.to_json()
                for citation in network.citations[offset_start:offset_end]
            ]

        if author:
            result = self.manager.query_citations(author=author, as_dict_list=True)

        if len(result) == 0:
            result = self.manager.query_citations(as_dict_list=True)

        return result

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

        #: dictionary of {tuple node: int id}
        self.node_nid = {}

        #: dictionary of {int id: tuple node}
        self.nid_node = {}

        #: dictionary of {str BEL: int id}
        self.bel_id = {}
        self.id_bel = {}

        #: The complete graph of all knowledge stored in the cache
        self.universe = BELGraph()

        self.universe_pmids = set()
        self.universe_authors = set()

        #: A dictionary from {int id: {tuple node: int degree}}
        self.node_degrees = {}

        #: A dictionary from {int network_id: {int network_id: float tanimoto node overlap}}
        self.overlap_cache = defaultdict(dict)

        self.eid_edge = {}
        self.edge_tuple_eid = {}

        if autocache:
            self.cache_networks()

        log.info('initialized dictionary service')

    def clear(self):
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

            if node in self.node_nid:
                continue

            node_id = data[ID]

            self.node_nid[node] = node_id
            self.nid_node[node_id] = node

            bel = node_to_bel(graph, node)
            self.id_bel[node_id] = bel
            self.bel_id[bel] = node_id

    def _relabel_notes_to_identifiers_helper(self, graph):
        if 'PYBEL_RELABELED' in graph.graph:
            log.warning('%s has already been relabeled', graph.name)
            return graph

        mapping = {node: self.node_nid[node] for node in graph.nodes_iter()}
        nx.relabel.relabel_nodes(graph, mapping, copy=False)

        graph.graph['PYBEL_RELABELED'] = True

        return graph

    def relabel_nodes_to_identifiers(self, graph, copy=True):
        """Relabels all nodes by their identifiers, in place. This function is a thin wrapper around
        :func:`networkx.relabel.relabel_nodes` with the module level variable :data:`node_nid` used as the mapping.

        :param pybel.BELGraph graph: A BEL Graph
        :param bool copy: Copy the graph first?
        :rtype: pybel.BELGraph
        """
        return self._relabel_notes_to_identifiers_helper(graph.copy() if copy else graph)

    def _add_network(self, network_id, force_reload=False, eager=False, maintain_universe=True):
        """Adds a network to the module-level cache from the underlying database

        :param int network_id: The identifier for the graph
        :param bool force_reload: Should the graphs be reloaded even if it has already been cached?
        :param bool eager: Should data be calculated/loaded eagerly?
        """
        if network_id in self.networks and not force_reload:
            log.info('tried re-adding graph [%s]', network_id)
            return self.networks[network_id]

        network = self.manager.get_network_by_id(network_id)

        try:
            log.debug('getting bytes from [%s]', network.id)
            graph = network.as_bel()
        except:
            log.exception("couldn't load from bytes [%s]", network.id)
            return

        t = time.time()

        log.info(
            'caching network [%s] %s v%s',
            network_id,
            graph.name,
            graph.version,
        )

        log.debug('parsing authors')
        parse_authors(graph)

        log.debug('inferring central dogma')
        infer_central_dogma(graph)

        log.debug('adding canonical names')
        add_canonical_names(graph)

        log.debug('updating graph node/edge indexes')
        self._update_indexes(graph)

        log.debug('calculating node degrees')
        self.node_degrees[network_id] = Counter(graph.degree())

        log.debug('caching PubMed identifiers')
        self.universe_pmids |= get_pubmed_identifiers(graph)

        if eager:
            log.debug('enriching citations')
            t = time.time()
            enrich_pubmed_citations(graph, manager=self.manager)
            log.debug('done enriching citations in %.2f seconds', time.time() - t)

        authors = get_authors(graph)
        log.debug('caching %d authors', len(authors))
        self.universe_authors |= authors

        if maintain_universe:
            log.debug('adding to the universe')
            left_full_join(self.universe, graph)

        self.networks[network_id] = graph

        log.info(
            'cached (%d nodes, %d edges) in %.2f seconds',
            graph.number_of_nodes(),
            graph.number_of_edges(),
            time.time() - t
        )

        return graph

    def cache_networks(self, force_reload=False, eager=False, maintain_universe=True):
        """This function needs to get all networks from the graph cache manager and make a dictionary

        :param bool force_reload: Should all graphs be reloaded even if they have already been cached?
        :param bool eager: Should difficult to preload features be calculated?
        """
        query = self.manager.session.query(Network).group_by(Network.name)

        if not force_reload and self.networks:
            query = query.filter(Network.id.notin_(self.networks))

        for network in query.having(func.max(Network.created)).order_by(Network.created.desc()).all():
            self._add_network(
                network.id,
                force_reload=force_reload,
                eager=eager,
                maintain_universe=maintain_universe
            )

        if maintain_universe:
            log.info(
                'universe has (%s nodes, %s edges)',
                self.universe.number_of_nodes(),
                self.universe.number_of_edges()
            )

    # Graph selection functions

    def get_graph_by_id(self, network_id=None):
        """Gets a network by its ID or super network if identifier is not specified

        :param int network_id: The internal ID of the network to get
        :return: A BEL Graph
        :rtype: pybel.BELGraph
        """
        if network_id is None:
            log.debug(
                'got universe (%s nodes, %s edges)',
                self.universe.number_of_nodes(),
                self.universe.number_of_edges()
            )
            return self.universe

        if network_id not in self.networks:
            log.debug('caching network [%s]', network_id)
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

    def get_node_hash(self, node):
        """Gets the hashes of the PyBEL node tuple

        :param tuple node: A PyBEL node tuple
        """
        return self.node_nid[node]

    def get_node_hashes(self, node_tuples):
        """Converts a list of BEL nodes to their node identifiers

        :param list[tuple] node_tuples: A list of PyBEL node tuples
        :rtype: list
        """
        return [
            self.get_node_hash(node)
            for node in node_tuples
        ]

    def get_node_tuple_by_hash(self, node_hash):
        """Returns the node tuple based on the node id

        :param node_hash: The node's identifier
        :return: A PyBEL node tuple
        :rtype: tuple
        """
        return self.nid_node.get(node_hash)

    def get_nodes_by_hashes(self, node_hashes):
        """Gets a list of node tuples from a list of ids

        :param list node_hashes: A list of node identifiers
        :rtype: list[tuple]
        """
        return [
            self.get_node_tuple_by_hash(node_hash)
            for node_hash in node_hashes
        ]

    def paths_tuples_to_ids(self, paths):
        """List of lists of node tuples

        :param list[list[tuple]] paths: list of lists (tuples)
        :rtype: list[list[tuple]]: list of lists (ids)
        """
        return [
            self.get_node_hashes(path)
            for path in paths
        ]

    def get_nodes_containing_keyword(self, keyword):
        """Gets a list with all cnames that contain a certain keyword adding to the duplicates their function

        :param str keyword: Search for nodes whose cnames have this as a substring
        :rtype: list[dict]
        """
        return [
            {"text": bel, "id": str(nid)}
            for bel, nid in self.bel_id.items()
            if keyword.lower() in bel.lower()
        ]

    def get_pubmed_containing_keyword(self, keyword):
        """Gets a list of PubMed identifiers that contain a certain keyword

        :param str keyword: Search for PubMed identifiers who have this as a substring
        :rtype: list[str]
        """
        return list(get_pmid_by_keyword(keyword, pubmed_identifiers=self.universe_pmids))

    def get_authors_containing_keyword(self, keyword):
        """Gets a list with authors that contain a certain keyword

        :param str keyword: Search for authors whose names have this as a substring
        :rtype: list[str]
        """
        # TODO switch to doing database lookup
        return get_authors_by_keyword(keyword, authors=self.universe_authors)

    def get_cname_by_node_hash(self, node_hash):
        """Gets the canonical name of a node

        :param node_hash: A BEL node identifier
        :rtype: str
        """
        node = self.get_node_tuple_by_hash(node_hash)
        return self.get_cname(node)

    def get_cname(self, node):
        """Gets the canonical name of a node

        :param tuple node: A BEL node
        :rtype: str
        """
        if CNAME in self.universe.node[node]:
            return self.universe.node[node][CNAME]

        self.universe.node[node][CNAME] = calculate_canonical_name(self.universe, node)
        return self.universe.node[node][CNAME]

    def get_top_degree(self, network_id, count=20):
        """Gets the nodes with the highest degrees

        :param int network_id: The network database identifier
        :param int count: The number of top degree nodes to get
        :rtype: dict[str,int]
        """
        if network_id not in self.node_degrees:
            log.info('lazy loading degrees for [%s]', network_id)
            graph = self.get_graph_by_id(network_id)
            self.node_degrees[network_id] = Counter(graph.degree())

        return {
            self.get_cname(node): v
            for node, v in self.node_degrees[network_id].most_common(count)
        }

    @lru_cache(maxsize=32)
    def get_top_pathologies(self, network_id, count=20):
        """Gets the top most frequent pathologies mentioned in a graph

        :param int network_id: The network database identifier
        :param int count: The number of most frequently mentioned pathologies to get
        :rtype: dict[str,int]
        """
        network = self.manager.get_network_by_id(network_id)
        graph = network.as_bel()
        cm = count_pathologies(graph)

        return {
            self.get_cname(node): v
            for node, v in cm.most_common(count)
        }

    @lru_cache(maxsize=32)
    def get_tree_annotations(self, graph):
        """Gets tree annotations for the given graph

        :param pybel.BELGraph graph: A BEL Graph
        :return: Annotations for the given graph
        :rtype: list[dict]
        """
        return get_tree_annotations(graph)

    def get_node_overlap(self, network_id):
        """Calculates overlaps to all other networks in the database

        :param int network_id: The network database identifier
        :return: A dictionary from {int network_id: float similarity} for this network to all other networks
        :rtype: collections.Counter[int, float]
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
        nodes = set(graph.nodes_iter())

        for other_network_id in other_network_ids:
            if other_network_id == network_id:
                continue

            if other_network_id in self.overlap_cache[network_id]:
                continue

            other_network = self.get_graph_by_id(other_network_id)
            other_network_nodes = set(other_network.nodes_iter())

            overlap = min_tanimoto_set_similarity(nodes, other_network_nodes)

            self.overlap_cache[network_id][other_network_id] = overlap
            self.overlap_cache[other_network_id][network_id] = overlap

        log.debug('Cached node overlaps for %s in %.2f seconds', graph, time.time() - t)

        return Counter(self.overlap_cache[network_id])

    def get_annotations_containing_keyword(self, keyword):
        """Gets annotation/value pairs for values for whom the search string is a substring

        :param str keyword: Search for annotations whose values have this as a substring
        :rtype: list[dict[str,str]
        """
        return get_annotations_containing_keyword(self.universe, keyword)

    def forget_network(self, network_id):
        """Removes all cached data from the given network id"""
        if network_id in self.networks:
            del self.networks[network_id]

        if network_id in self.node_degrees:
            del self.node_degrees[network_id]

        if network_id in self.overlap_cache:
            del self.overlap_cache[network_id]

            for k in self.overlap_cache:
                if network_id in self.overlap_cache[k]:
                    del self.overlap_cache[k][network_id]
