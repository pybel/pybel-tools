# -*- coding: utf-8 -*-

import logging
import warnings

import obonet

from .node_annotator import NodeAnnotator

__all__ = [
    'GOAnnotator',
]

log = logging.getLogger(__name__)

url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'


class GOAnnotator(NodeAnnotator):
    """Annotates GO entries"""

    def __init__(self, preload=True):
        """
        :param bool preload: Should the data be preloaded?
        """
        super(GOAnnotator, self).__init__([
            'GOBP',
            'GOBPID',
            'GOCC',
            'GOCCID',
            'GOMF',
            'GOMFID',
        ])

        warnings.warn("deprecated. Use bio2bel_go.GOAnnotator instead", DeprecationWarning)

        #: A dictionary of {str go term/id: str description}
        self.descriptions = {}

        if preload:
            self.download_successful = self.download()

    # OVERRIDES
    def get_description(self, name):
        return self.descriptions.get(name)

    def download(self):
        """Downloads the OBO file for Gene Ontology. Returns true on success."""
        try:
            self.graph = obonet.read_obo(url)
        except:
            log.exception('Unable to download Gene Ontology OBO')
            return False

        self.id_to_name = {id_: data['name'] for id_, data in self.graph.nodes(data=True)}
        self.name_to_id = {data['name']: id_ for id_, data in self.graph.nodes(data=True)}

        for id_, data in self.graph.nodes(data=True):
            self.descriptions[data['name']] = data['def']
            self.descriptions[id_] = data['def']

        return True
