# -*- coding: utf-8 -*-

import json
import logging
from collections import Iterable

import numpy as np

from pybel.struct import union
from pybel.utils import list2tuple
from .pipeline import Pipeline
from .selection import get_subgraph
from .selection.induce_subgraph import (
    NONNODE_SEED_TYPES, SEED_TYPE_ANNOTATION, SEED_TYPE_INDUCTION, SEED_TYPE_NEIGHBORS, SEED_TYPE_SAMPLE,
)

log = logging.getLogger(__name__)

SEED_TYPE_KEY = 'type'
SEED_DATA_KEY = 'data'


class Query:
    """Wraps a query over the network store"""

    def __init__(self, network_ids=None, seeding=None, pipeline=None):
        """
        :param iter[int] network_ids: Database network identifiers identifiers
        :param Optional[list[dict]] seeding:
        :param Optional[Pipeline] pipeline: Instance of a pipeline
        """
        if not network_ids:
            self.network_ids = []
        else:
            if not isinstance(network_ids, Iterable):
                raise TypeError('network identifiers is not list: {}'.format(network_ids))

            network_ids = list(network_ids)

            if any(not isinstance(entry, int) for entry in network_ids):
                raise TypeError('network identifiers entry is not int: {}'.format(network_ids))

            self.network_ids = network_ids

        self.seeding = seeding or []
        self.pipeline = pipeline or Pipeline()

    def append_network(self, network_id):
        """Adds a network to this query

        :param int network_id: The database identifier of the network
        """
        self.network_ids.append(network_id)

    def append_seed(self, seed_type, data):
        """Add a seeding

        :param str seed_type:
        :param data:
        """
        self.seeding.append({
            SEED_TYPE_KEY: seed_type,
            SEED_DATA_KEY: data
        })

    def append_seeding_induction(self, data):
        """Adds a seed induction method

        :param list[tuple] data: A list of PyBEL node tuples
        """
        self.append_seed(SEED_TYPE_INDUCTION, data)

    def append_seeding_neighbors(self, data):
        """Adds a seed by neighbors

        :param list[tuple] data:
        """
        self.append_seed(SEED_TYPE_NEIGHBORS, data)

    def append_seeding_annotation(self, annotation, values):
        """Adds a seed induction method for single annotation's values

        :param str annotation: The annotation to filter by
        :param set[str] values: The values of the annotation to keep
        """
        self.append_seed(SEED_TYPE_ANNOTATION, {
            'annotations': {
                annotation: values
            }
        })

    def append_seeding_sample(self):
        """Adds seed induction methods.

        Kwargs can have ``number_edges`` or ``number_seed_nodes``.
        """
        self.append_seed(SEED_TYPE_SAMPLE, {
            'seed': np.random.randint(0, np.iinfo('i').max)
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

        if not self.seeding:
            return self.pipeline.run(query_universe, universe=query_universe)

        subgraphs = []
        for seed in self.seeding:
            subgraph = get_subgraph(
                query_universe,
                seed_method=seed[SEED_TYPE_KEY],
                seed_data=seed[SEED_DATA_KEY]
            )

            if subgraph is None:
                log.debug('Seed returned empty graph: %s', seed)
                continue

            subgraphs.append(subgraph)

        if not subgraphs:
            return

        graph = union(subgraphs)

        return self.pipeline.run(graph, universe=query_universe)

    def seeding_to_jsons(self):
        """Returns seeding JSON as a string

        :rtype: str
        """
        return json.dumps(self.seeding)

    def to_json(self):
        """Returns this query as a JSON object

        :rtype: dict
        """
        rv = {
            'network_ids': list(self.network_ids),
        }

        if self.seeding:
            rv['seeding'] = self.seeding

        if self.pipeline:
            rv['pipeline'] = self.pipeline.to_json()

        return rv

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
        rv = Query(network_ids=d['network_ids'])

        if 'seeding' in d:
            rv.seeding = process_seeding(d['seeding'])

        if 'pipeline' in d:
            rv.pipeline.protocol = d['pipeline']

        return rv


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
