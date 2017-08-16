# -*- coding: utf-8 -*-

from __future__ import print_function

from pybel.constants import HAS_PRODUCT, HAS_REACTANT, HAS_VARIANT, HAS_COMPONENT, TRANSCRIBED_TO, TRANSLATED_TO

IS_PRODUCT_OF = 'isProductOf'
IS_REACTANT_OF = 'isReactantOf'
IS_VARIANT_OF = 'isVariantOf'
IS_COMPONENT_OF = 'isComponentOf'
TRANSCRIBED_FROM = 'transcribedFrom'
TRANSLATED_FROM = 'translatedFrom'

INFERRED_INVERSE = {
    HAS_PRODUCT: IS_PRODUCT_OF,
    HAS_REACTANT: IS_REACTANT_OF,
    HAS_VARIANT: IS_VARIANT_OF,
    HAS_COMPONENT: IS_COMPONENT_OF,
    TRANSCRIBED_TO: TRANSCRIBED_FROM,
    TRANSLATED_TO: TRANSLATED_FROM
}

abstract_url_fmt = "http://togows.dbcls.jp/entry/ncbi-pubmed/{}/abstract"
title_url_fmt = "http://togows.dbcls.jp/entry/ncbi-pubmed/{}/title"
#: SO gives short citation information
so_url_fmt = "http://togows.dbcls.jp/entry/ncbi-pubmed/{}/so"
citation_format = 'SET Citation = {{"PubMed","{}","{}"}}'
evidence_format = 'SET Evidence = "{}"'


def pubmed(name, identifier):
    return citation_format.format(name.replace('\n', ''), identifier)


def print_set_pubmed(pmid, title=None, file=None):
    if title:
        print(pubmed(title, pmid), file=file)
    else:
        print('SET Citation = {{"PubMed","{}"}}'.format(pmid), file=file)


def print_set_evidence(evidence, file=None):
    print(evidence_format.format(evidence), file=file)


def print_set_confidence(confidence, file=None):
    print('SET Confidence = "{}"'.format(confidence), file=file)


CNAME = 'cname'
PUBMED = 'PubMed'
DATA_WEIGHT = 'weight'

# Resources

#: URL For HGNC Gene Families Memberships BEL document
GENE_FAMILIES = 'https://arty.scai.fraunhofer.de/artifactory/bel/knowledge/hgnc-gene-family-membership/hgnc-gene-family-membership-20170710.bel'
NAMED_COMPLEXES = 'http://resources.openbel.org/belframework/20150611/resource/named-complexes.bel'

#: Points to the env variable name for PyBEL resources
PYBEL_RESOURCES_ENV = 'PYBEL_RESOURCES_BASE'

#: Points to the env variable for ownCloud resources
OWNCLOUD_ENV = 'OWNCLOUD_BASE'

#: Points to the env variable for the biological model store repository
BMS_BASE = 'BMS_BASE'

DEFAULT_SERVICE_URL = 'https://pybel.scai.fraunhofer.de'
