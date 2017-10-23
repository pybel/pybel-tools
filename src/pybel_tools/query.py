# -*- coding: utf-8 -*-

import json
import logging

from pybel.struct import union
from pybel.utils import list2tuple
from .pipeline import Pipeline
from .selection import get_subgraph
from .selection.induce_subgraph import (
    NONNODE_SEED_TYPES, SEED_TYPE_ANNOTATION, SEED_TYPE_INDUCTION,
    SEED_TYPE_NEIGHBORS,
)

log = logging.getLogger(__name__)

SEED_TYPE_KEY = 'type'
SEED_DATA_KEY = 'data'


class Query:
    """Wraps a query over the network store"""

    def __init__(self, network_ids, seed_list=None, pipeline=None):
        """
        :param int or list[int] or set[int] or tuple[int] network_ids:
        :param list[dict] seed_list:
        :param Pipeline pipeline: Instance of a pipeline
        """
        self.network_ids = []
        self.seeds = []

        if isinstance(network_ids, int):
            self.append_network(network_ids)
        elif isinstance(network_ids, (list, set, tuple)):
            self.network_ids.extend(int(network_id) for network_id in network_ids)
        else:
            raise TypeError(network_ids)

        if seed_list:
            self.seeds.extend(seed_list)

        self.pipeline = pipeline if pipeline is not None else Pipeline()

    def append_network(self, network_id):
        """Adds a network to this query

        :param int network_id: The database identifier of the network
        """
        self.network_ids.append(network_id)

    def add_seed(self, seed_type, data):
        self.seeds.append({
            SEED_TYPE_KEY: seed_type,
            SEED_DATA_KEY: data
        })

    def add_seed_induction(self, data):
        """Adds a seed induction method

        :param list[tuple] data: A list of PyBEL node tuples
        """
        self.add_seed(SEED_TYPE_INDUCTION, data)

    def add_seed_neighbors(self, data):
        """Adds a seed by neighbors

        :param list[tuple] data:
        """
        self.add_seed(SEED_TYPE_NEIGHBORS, data)

    def add_seed_annotation(self, annotation, values):
        """Adds a seed induction method for single annotation's values

        :param str annotation: The annotation to filter by
        :param set[str] values: The values of the annotation to keep
        """
        self.add_seed(SEED_TYPE_ANNOTATION, {
            'annotations': {
                annotation: values
            }
        })

    def append_pipeline(self, name, *args, **kwargs):
        """Adds an entry to the pipeline

        :param str name: The name of the function
        :return: This pipeline for fluid query building
        :rtype: Pipeline
        """
        return self.pipeline.append(name, *args, **kwargs)

    def run(self, manager):
        """Runs this query and returns the resulting BEL graph

        :param pybel.manager.Manager manager: A cache manager
        :rtype: Optional[pybel.BELGraph]
        """
        log.debug('query universe consists of networks: %s', self.network_ids)

        if not self.network_ids:
            return

        query_universe = manager.get_graph_by_ids(self.network_ids)

        log.debug(
            'query universe has %d nodes/%d edges',
            query_universe.number_of_nodes(),
            query_universe.number_of_edges()
        )

        # parse seeding stuff

        if not self.seeds:
            return self.pipeline.run(query_universe, universe=query_universe)

        subgraphs = []
        for seed in self.seeds:
            subgraph = get_subgraph(
                query_universe,
                seed_method=seed[SEED_TYPE_KEY],
                seed_data=seed[SEED_DATA_KEY]
            )

            if subgraph is None:
                log.debug('Seed returned empty graph: %s', seed)
                continue

            # TODO streamline this logging... maybe put in get_subgraph function
            log.debug(
                'Subgraph coming from %s (seed type) %s (data) contains %d nodes and %d edges',
                seed[SEED_DATA_KEY],
                seed[SEED_TYPE_KEY],
                subgraph.number_of_nodes(),
                subgraph.number_of_edges()
            )

            subgraphs.append(subgraph)

        if not subgraphs:
            return

        graph = union(subgraphs)

        log.debug(
            'Number of nodes/edges in query before running pipeline: %d nodes, %d edges )',
            graph.number_of_nodes(),
            graph.number_of_edges(),
        )

        return self.pipeline.run(graph, universe=query_universe)

    def seeding_to_json(self):
        """Returns seeding JSON

        :rtype: list[dict]
        """
        return [
            {
                SEED_TYPE_KEY: seed[SEED_TYPE_KEY],
                SEED_DATA_KEY: seed[SEED_DATA_KEY]
            }
            for seed in self.seeds
        ]

    def seeding_to_jsons(self):
        """Returns seeding JSON as a string

        :rtype: str
        """
        return json.dumps(self.seeding_to_json())

    def to_json(self):
        """Returns this query as a JSON object

        :rtype: dict
        """
        return {
            'network_ids': list(self.network_ids),
            'seeding': self.seeding_to_json(),
            'pipeline': self.pipeline.to_json()
        }

    def to_jsons(self):
        """Returns this query as a stringified JSON object

        :rtype: str
        """
        return json.dumps(self.to_json())

    @staticmethod
    def from_jsons(s):
        """Loads a query from a stringified JSON object

        :param str s: A stringified JSON query
        :rtype: Query
        """
        return Query.from_json(json.loads(s))

    @staticmethod
    def from_json(d):
        """Loads a query from a JSON dictionary

        :param dict d: A JSON dictionary
        :rtype: Query
        """
        return Query(
            network_ids=d['network_ids'],
            seed_list=process_seeding(d['seeding']),
            pipeline=Pipeline.from_json(d['pipeline']),
        )


def process_seeding(seeds):
    """Makes sure nodes are tuples and not lists once back in"""
    return [
        {
            SEED_TYPE_KEY: seed[SEED_TYPE_KEY],
            SEED_DATA_KEY: [
                list2tuple(node)
                for node in seed[SEED_DATA_KEY]
            ]
            if seed[SEED_TYPE_KEY] not in NONNODE_SEED_TYPES else seed[SEED_DATA_KEY]
        }
        for seed in seeds
    ]
