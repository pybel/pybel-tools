# -*- coding: utf-8 -*-

from __future__ import print_function

from pybel.constants import HAS_COMPONENT, HAS_PRODUCT, HAS_REACTANT, HAS_VARIANT, TRANSCRIBED_TO, TRANSLATED_TO
from pybel.resources.constants import citation_format, evidence_format

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

CNAME = 'cname'
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
