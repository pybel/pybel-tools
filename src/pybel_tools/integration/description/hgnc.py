# -*- coding: utf-8 -*-

import logging
import time
import warnings

import pandas as pd

from .node_annotator import NodeAnnotator
from ...document_utils import get_entrez_gene_data
from ...summary import get_names_by_namespace
from ...utils import grouper

__all__ = [
    'HGNCAnnotator',
]

log = logging.getLogger(__name__)

#: Looks up the HGNC Gene Symbol to Entrez Gene Identifier mapping from HGNC
HGNC_ENTREZ_URL = 'http://www.genenames.org/cgi-bin/download?col=gd_app_sym&col=gd_pub_eg_id&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit'


class HGNCAnnotator(NodeAnnotator):
    """Annotates the labels and descriptions of Genes with HGNC identifiers using a mapping provided by HGNC
    and then the Entrez Gene Service.
    """

    def __init__(self, preload=True):
        """
        :param bool preload: Should the data be pre-downloaded?
        """
        super(HGNCAnnotator, self).__init__('HGNC')

        warnings.warn("deprecated. Use bio2bel_hgnc.HGNCAnnotator instead", DeprecationWarning)

        #: A dictionary of {str hgnc gene symbol: str description}
        self.descriptions = {}
        #: A dictionary of {str hgnc gene symbol: str label}
        self.labels = {}

        if preload:
            self.download_successful = self.load_hgnc_entrez_map()

    # OVERRIDES
    def populate_by_graph(self, graph):
        """Downloads the gene information only for genes in the given graph

        :param pybel.BELGraph graph: A BEL graph
        """
        hgnc_symbols = list(get_names_by_namespace(graph, self.namespace))
        self.populate_constrained(hgnc_symbols)

    # OVERRIDES
    def get_description(self, name):
        return self.descriptions.get(name)

    # OVERRIDES
    def get_label(self, name):
        return self.labels.get(name)

    def load_hgnc_entrez_map(self):
        """Preloads the HGNC-Entrez map"""
        log.info('Downloading HGNC-Entrez map from %s', HGNC_ENTREZ_URL)

        try:
            df = pd.read_csv(HGNC_ENTREZ_URL, sep='\t')
        except:
            log.exception('Unable to download HGNC Entrez map')
            return False

        self.hgnc_entrez = {
            hgnc: entrez_id
            for _, hgnc, entrez_id in df[['Approved Symbol', 'Entrez Gene ID']].itertuples()
        }

        self.entrez_hgnc = {v: k for k, v in self.hgnc_entrez.items()}

        return True

    def map_entrez_ids(self, entrez_ids):
        """Maps a list of Entrez Gene Identifiers to HGNC Gene Symbols"""
        return [self.entrez_hgnc[entrez] for entrez in entrez_ids]

    def map_hgnc(self, hgnc_symbols):
        """Maps a list of HGNC Gene Symbols to Entrez Gene Identifiers"""
        return [self.hgnc_entrez[hgnc] for hgnc in hgnc_symbols if hgnc in self.hgnc_entrez]

    def get_unpopulated_entrez(self, entrez_ids):
        """Gets the Entrez Gene Identifiers from this list that aren't already cached"""
        hgnc_symbols = self.map_entrez_ids(entrez_ids)
        missing_hgnc_symbols = set(hgnc_symbols) - set(self.descriptions)
        return self.map_hgnc(missing_hgnc_symbols)

    def populate(self, entrez_ids, group_size=200, sleep_time=1):
        """Download the descriptions from Entrez Gene Service for a given list of Entrez Gene Identifiers

        :param iter entrez_ids: An iterable of Entrez Gene Identifiers
        :param int group_size: The number of entrez gene id's to send per query
        :param int sleep_time: The number of seconds to sleep between queries
        """
        unpopulated_entrez = self.get_unpopulated_entrez(entrez_ids)

        for entrez_ids_group in grouper(group_size, unpopulated_entrez):
            for entrez_id, data in get_entrez_gene_data(entrez_ids_group).items():
                try:
                    hgnc = self.entrez_hgnc[int(entrez_id)]

                    self.labels[hgnc] = data['description']
                    self.descriptions[hgnc] = data['summary']
                except Exception:
                    log.warning('Missing for EGID: %s. Dat: %s', entrez_id, data)

            time.sleep(sleep_time)

    def populate_unconstrained(self, group_size=200, sleep_time=1):
        """Downloads all descriptions for all Entrez Gene Identifiers"""

        #: This variable keeps track of when the data was downloaded
        self.download_time = time.asctime()

        self.populate(
            entrez_ids=self.hgnc_entrez.values(),
            group_size=group_size,
            sleep_time=sleep_time,
        )

    def populate_constrained(self, hgnc_symbols, group_size=200, sleep_time=1):
        """Downloads the gene information only for genes in the list of HGNC Gene Symbols"""
        entrez_ids = self.map_hgnc(hgnc_symbols)

        self.populate(
            entrez_ids=entrez_ids,
            group_size=group_size,
            sleep_time=sleep_time,
        )
