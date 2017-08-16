# -*- coding: utf-8 -*-

"""

Chemical Similarity Enrichment
------------------------------
1. Extract all chemicals from ChEBI from BEL network

2. Get database data
    a. Get ChEBI names to ChEBI identifiers: ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz
    b. Get ChEBI identifiers to InChI: ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv
3. Map all ChEBI to InChI
4. Calculate all pairwise relationships between stuff in network and add associations for certain threshold.

Seeding by Chemical Similarity
------------------------------
- Given an InChI string, calculate all similarities to all compounds annotated in the network. Seed on node list
  of similar ones
- Calculate enrichment of similar compounds effects on certain networks (using derived NeuroMMSig algorithm)

"""

import xml.etree.ElementTree as ET

import requests

from ..summary.node_summary import get_names_by_namespace

#: This key gets put in the node data dictionary to show that the node has been annotated with an inchi key
INCHI_KEY = 'inchikey'

CHEBI_NAME_ENDPOINT = 'http://www.ebi.ac.uk/webservices/chebi/2.0/test/getLiteEntity?searchCategory=CHEBI+NAME&maximumResults=200&starsCategory=ALL&search='
CHEBI_INFO_ENDPOINT = 'http://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiid='
CHEBI_ID_XPATH = '{http://www.ebi.ac.uk/webservices/chebi}chebiId'
CHEBI_NAME_XPATH = '{http://www.ebi.ac.uk/webservices/chebi}chebiAsciiName'

PYBEL_ID_KEY = 'pybel_id'
PYBEL_INCHI_KEY = 'pybel_inchi'


def annotate_chebi_names(graph):
    """Annotates InChI keys to all nodes with the ChEBI namespace
    
    :param pybel.BELGraph graph: A BEL Graph 
    """
    names = get_names_by_namespace(graph, 'chebi')

    name_id = {}
    id_inchi = {}

    for name in names:
        res = requests.get(CHEBI_NAME_ENDPOINT + name)
        tree = ET.fromstring(res.content)
        data = {x.tag: x.text for x in tree[0][0][0][0]}

        name_id[data[CHEBI_NAME_XPATH]] = data[CHEBI_ID_XPATH]

        # TODO
        # for chebiid in set(name_id.values()):
        #    res_info = requests.get(CHEBI_INFO_ENDPOINT + chebiid)


def annotate_chebi_ids(graph):
    """Annotates InChI keys to all nodes with the ChEBI-ID namespace

    :param pybel.BELGraph graph: A BEL Graph 
    """
    raise NotImplementedError


def annotate_chembl_ids(graph):
    """Annotates InChI keys to all nodes with the ChEMBL namespace

    :param pybel.BELGraph graph: A BEL Graph 
    """
    raise NotImplementedError


def annotate_inchi(graph):
    """Annotates InChI keys to all nodes identified with InChI
    
    :param pybel.BELGraph graph: A BEL Graph 
    """
    raise NotImplementedError
