# -*- coding: utf-8 -*-

"""This file contains the curated resource dictionary"""

import logging
import os

from artifactory import ArtifactoryPath

from .definition_utils import get_bel_resource_hash
from .utils import get_iso_8601_date

log = logging.getLogger(__name__)

ARTY_BASE = 'https://arty.scai.fraunhofer.de/artifactory/bel/'
ARTY_NS = ARTY_BASE + 'namespace/'
ARTY_ANNO = ARTY_BASE + 'annotation/'
ARTY_BEL = ARTY_BASE + 'knowledge/'


def get_arty_namespace_module(namespace):
    return '{0}{1}/'.format(ARTY_NS, namespace)


def get_arty_namespace(namespace, version):
    return '{}-{}.belns'.format(namespace, version)


def get_arty_namespace_url(namespace, version):
    """Gets a BEL namespace file from artifactory given the name and version"""
    return '{}{}'.format(get_arty_namespace_module(namespace), get_arty_namespace(namespace, version))


def get_arty_annotation_module(module_name):
    return '{0}{1}/'.format(ARTY_ANNO, module_name)


def get_arty_annotation(module_name, version):
    return '{}-{}.belanno'.format(module_name, version)


def get_arty_annotation_url(module_name, version):
    """Gets a BEL annotation file from artifactory given the name and version"""
    return '{}{}'.format(get_arty_annotation_module(module_name), get_arty_annotation(module_name, version))


def get_arty_knowledge_module(module_name):
    return '{0}{1}/'.format(ARTY_BEL, module_name)


def get_arty_knowledge(module_name, version):
    return '{}-{}.bel'.format(module_name, version)


def get_arty_knowledge_url(module_name, version):
    """Gets a BEL knowledge file from artifactory given the name and version"""
    return '{}{}'.format(get_arty_knowledge_module(module_name), get_arty_knowledge(module_name, version))


def get_today_arty_namespace(module_name):
    """Gets the right name for the next version of the namespace"""
    return get_arty_namespace(module_name, get_iso_8601_date())


def get_today_arty_annotation(module_name):
    """Gets the right name for the next version of the annotation"""
    return get_arty_annotation(module_name, get_iso_8601_date())


def get_today_arty_knowledge(knowledge):
    return get_arty_knowledge(knowledge, get_iso_8601_date())


def get_arty_auth():
    """Gets the arty authentication tuple from the environment variables ``ARTY_USERNAME`` and ``ARTY_PASSWORD``,
    respectively.

    :rtype: tuple[str]
    """
    return os.environ['ARTY_USERNAME'], os.environ['ARTY_PASSWORD']


def _get_path_helper(module_name, getter):
    """Helps get the Artifactory path for a certain module

    :param str module_name: The name of the module
    :param types.FunctionType getter: The function that gets the modules from the Artifactory repository
    :rtype: artifactory.ArtifactoryPath
    """
    return ArtifactoryPath(getter(module_name))


def get_namespace_history(module_name):
    """Gets the Artifactory path for a namespace module

    :param str module_name: The name of the namespace module
    :rtype: artifactory.ArtifactoryPath
    """
    return _get_path_helper(module_name, get_arty_namespace_module)


def get_annotation_history(module_name):
    """Gets the Artifactory path for an annotation module

    :param str module_name: The name of the annotation module
    :rtype: artifactory.ArtifactoryPath
    """
    return _get_path_helper(module_name, get_arty_annotation_module)


def get_knowledge_history(module_name):
    """Gets the Artifactory path for a knowledge module

    :param str module_name: The name of the knowledge module
    :rtype: artifactory.ArtifactoryPath
    """
    return _get_path_helper(module_name, get_arty_knowledge_module)


def _get_latest_arty_helper(module_name, getter):
    """Helps get the latest path for a given BEL module by paremetrizing the getter"""
    path = _get_path_helper(module_name, getter)
    mp = max(path)
    return mp.as_posix()


def get_latest_arty_namespace(module_name):
    """Gets the latest path for this BEL namespace module. For historical reasons, some of these are not the same
    as the keyword. For example, the module name for HGNC is ``hgnc-human-genes`` due to the Selventa nomenclature.
    See https://arty.scai.fraunhofer.de/artifactory/bel/namespace/ for the entire manifest of available namespaces.

    :param str module_name: The BEL namespace module name
    :return: The URL of the latest version of this namespace
    :rtype: str
    """
    return _get_latest_arty_helper(module_name, get_arty_namespace_module)


def get_latest_arty_annotation(module_name):
    """Gets the latest path for this BEL annotation module

    :param str module_name: The BEL annotation module name
    :return: The URL of the latest version of this annotation
    :rtype: str
    """
    return _get_latest_arty_helper(module_name, get_arty_annotation_module)


def get_latest_arty_knowledge(module_name):
    """Gets the latest path for this BEL annotation module

    :param str module_name: The BEL knowledge module name
    :return: The URL of the latest version of this knowledge document
    :rtype: str
    """
    return _get_latest_arty_helper(module_name, get_arty_knowledge_module)


def _deploy_helper(filename, module_name, get_module, get_today_fn, hash_check=True, auth=None):
    """Deploys a file to the Artifactory BEL namespace cache

    :param str filename: The physical path
    :param str module_name: The name of the module to deploy to
    :param tuple[str] auth: A pair of (str username, str password) to give to the auth keyword of the constructor of
                            :class:`artifactory.ArtifactoryPath`. Defaults to the result of :func:`get_arty_auth`.
    :return: The resource path, if it was deployed successfully, else none.
    :rtype: str
    """
    path = ArtifactoryPath(
        get_module(module_name),
        auth=get_arty_auth() if auth is None else auth
    )
    path.mkdir(exist_ok=True)

    if hash_check:
        deployed_semantic_hashes = {
            get_bel_resource_hash(subpath.as_posix())
            for subpath in path
        }

        semantic_hash = get_bel_resource_hash(filename)

        if semantic_hash in deployed_semantic_hashes:
            return  # Don't deploy if it's already uploaded

    target = path / get_today_fn(module_name)
    target.deploy_file(filename)

    log.info('deployed %s', module_name)

    return target.as_posix()


def deploy_namespace(filename, module_name, hash_check=True, auth=None):
    """Deploys a file to the Artifactory BEL namespace cache

    :param str filename: The physical path
    :param str module_name: The name of the module to deploy to
    :param bool hash_check: Ensure the hash is unique before deploying
    :param tuple[str] auth: A pair of (str username, str password) to give to the auth keyword of the constructor of
                            :class:`artifactory.ArtifactoryPath`. Defaults to the result of :func:`get_arty_auth`.
    :return: The resource path, if it was deployed successfully, else none.
    :rtype: str
    """
    return _deploy_helper(
        filename,
        module_name,
        get_arty_namespace_module,
        get_today_arty_namespace,
        hash_check=hash_check,
        auth=auth
    )


def deploy_knowledge(filename, module_name, auth=None):
    """Deploys a file to the Artifactory BEL knowledge cache

    :param str filename: The physical file path
    :param str module_name: The name of the module to deploy to
    :param tuple[str] auth: A pair of (str username, str password) to give to the auth keyword of the constructor of
                            :class:`artifactory.ArtifactoryPath`. Defaults to the result of :func:`get_arty_auth`.
    :return: The resource path, if it was deployed successfully, else none.
    :rtype: str
    """
    return _deploy_helper(
        filename,
        module_name,
        get_arty_knowledge_module,
        get_today_arty_knowledge,
        hash_check=False,
        auth=auth
    )


def deploy_annotation(filename, module_name, hash_check=True, auth=None):
    """Deploys a file to the Artifactory BEL annotation cache

    :param str filename: The physical file path
    :param str module_name: The name of the module to deploy to
    :param bool hash_check: Ensure the hash is unique before deploying
    :param tuple[str] auth: A pair of (str username, str password) to give to the auth keyword of the constructor of
                            :class:`artifactory.ArtifactoryPath`. Defaults to the result of :func:`get_arty_auth`.
    :return: The resource path, if it was deployed successfully, else none.
    :rtype: str
    """
    return _deploy_helper(
        filename,
        module_name,
        get_arty_annotation_module,
        get_today_arty_annotation,
        hash_check=hash_check,
        auth=auth
    )


HGNC_HUMAN_GENES = get_arty_namespace_url('hgnc-human-genes', '20170511')
CHEBI = get_arty_namespace_url('chebi', '20170511')
CHEBI_IDS = get_arty_namespace_url('chebi-ids', '20170511')
HGNC_GENE_FAMILIES = get_arty_namespace_url('hgnc-gene-families', '20170515')
CONFIDENCE = get_arty_annotation_url('confidence', '1.0.0')
MESHD = get_arty_annotation_url('mesh-diseases', '20170511')
NEUROMMSIG = get_arty_annotation_url('neurommsig', '1.0.1')
NIFT = get_arty_namespace_url('imaging-ontology', '1.0.0')

default_namespace_spec = [
    ('ADO', 'alzheimer-disease-ontology', '1.0.2'),
    ('AFFX', 'affx-probeset-ids', '20170511'),
    ('BRCO', 'brain-region-ontology', '1.0.0'),
    ('CHEBI', 'chebi', '20170511'),
    ('CHEBIID', 'chebi-ids', '20170511'),
    ('CTO', 'clinical-trial-ontology', '1.0.0'),
    ('DO', 'disease-ontology', '20170511'),
    ('EGID', 'entrez-gene-ids', '20170511'),
    ('EPT', 'epilepsy-terminology', '1.0.0'),
    ('FlyBase', 'flybase', '20170508'),
    ('GOBP', 'go-biological-process', '20170511'),
    ('GOCC', 'go-cellular-component', '20170511'),
    ('GFAM', 'hgnc-gene-families', '20170515'),
    ('HGNC', 'hgnc-human-genes', '20170511'),
    ('NIFT', 'imaging-ontology', '1.0.0'),
    ('NTN', 'nutrition', '1.0.0'),
    ('MESHCS', 'mesh-cell-structures', '20170511'),
    ('MESHD', 'mesh-diseases', '20170511'),
    ('MESHPP', 'mesh-processes', '20170511'),
    ('MGI', 'mgi-mouse-genes', '20170511'),
    ('PTS', 'neurodegeneration-pathways', '20170511'),
    ('PDO', 'parkinson-disease-ontology', '20170511'),
    ('RGD', 'rgd-rat-genes', '20170511'),
    ('SCOM', 'selventa-named-complexes', '20170511'),
    ('SFAM', 'selventa-protein-families', '20170511'),
    ('SP', 'swissprot', '20170511'),
]

default_namespaces = {
    keyword: get_arty_namespace_url(namespace, version)
    for keyword, namespace, version in default_namespace_spec
}

# See: https://gist.github.com/lsauer/1312860
DBSNP_PATTERN = 'rs[0-9]+'
EC_PATTERN = '(\d+|\-)\.((\d+)|(\-))\.(\d+|\-)(\.(n)?(\d+|\-))*'
INCHI_PATTERN = '^((InChI=)?[^J][0-9BCOHNSOPrIFla+\-\(\)\\\/,pqbtmsih]{6,})$'

default_namespace_patterns = {
    'dbSNP': DBSNP_PATTERN,
    'EC': EC_PATTERN,
    'InChI': INCHI_PATTERN
}
default_annotation_spec = [
    ('Anatomy', 'anatomy', '20170511'),
    ('Cell', 'cell', '20170511'),
    ('CellLine', 'cell-line', '20170511'),
    ('CellStructure', 'cell-structure', '20170511'),
    ('Confidence', 'confidence', '1.0.0'),
    ('Disease', 'disease', '20170511'),
    ('Gender', 'gender', '1.0.0'),
    ('MeSHAnatomy', 'mesh-anatomy', '20170511'),
    ('MeSHDisease', 'mesh-diseases', '20170511'),
    ('Subgraph', 'neurommsig', '1.0.1'),
    ('SNPO', 'snpo', '20170425'),
    ('Species', 'species-taxonomy-id', '20170511'),
    ('TextLocation', 'text-location', '1.0.0'),
]
belief_demo_prefix_1 = 'http://belief-demo.scai.fraunhofer.de/openbel/repository/namespaces/'
belief_demo_namespaces_1 = {
    'CHEMBL': 'chembl-names.belns',
    'CHEMBLID': 'chembl-ids.belns',
    'LMSD': 'LMSD.belns',
}

belief_demo_prefix_2 = 'http://belief-demo.scai.fraunhofer.de/BeliefDashboard/dicten/namespaces/'
belief_demo_namespaces_2 = {
    'PMIBP': 'pmibp.belns',
    'PMICHEM': 'pmichem.belns',
    'PMICOMP': 'pmicomp.belns',
    'PMIDIS': 'pmidis.belns',
    'PMIPFAM': 'pmipfam.belns',
}

default_annotations = {
    keyword: get_arty_annotation_url(annotation, version)
    for keyword, annotation, version in default_annotation_spec
}


def deploy_directory(directory, auth=None):
    """Uploads all stuff from a directory to artifactory

    :param str directory: the path to a directory
    :param tuple[str] auth: A pair of (str username, str password) to give to the auth keyword of the constructor of
                            :class:`artifactory.ArtifactoryPath`. Defaults to the result of :func:`get_arty_auth`.
    """
    for file in os.listdir(directory):
        full_path = os.path.join(directory, file)

        if file.endswith('.belanno'):
            name = file[:-8]
            log.info('Uploading annotation %s', full_path)
            deploy_annotation(full_path, name, auth=auth)
        elif file.endswith('.belns'):
            name = file[:-6]
            log.info('Uploading namespace %s', full_path)
            deploy_namespace(full_path, name, auth=auth)
        elif file.endswith('.bel'):
            name = file[:-4]
            log.info('Uploading knowledge %s', full_path)
            deploy_knowledge(full_path, name, auth=auth)
