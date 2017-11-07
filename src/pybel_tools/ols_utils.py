# -*- coding: utf-8 -*-

from collections import defaultdict

from ols_client import OlsClient
from pybel.constants import CITATION_TYPE_URL, IS_A, belns_encodings
from pybel.parser.language import rev_abundance_labels
from pybel.resources.arty import (
    get_latest_arty_namespace, get_today_arty_annotation, get_today_arty_knowledge,
    get_today_arty_namespace,
)
from pybel.resources.definitions import write_annotation, write_namespace
from pybel.resources.deploy import deploy_annotation, deploy_knowledge, deploy_namespace
from pybel.utils import ensure_quotes
from pybel_tools.document_utils import write_boilerplate

function_to_encoding = defaultdict(list)
for enc, functions in belns_encodings.items():
    for fn in functions:
        function_to_encoding[fn].append(enc)


class OlsNamespaceOntology:
    """Wraps the functions needed to use the OLS to generate and deploy BEL namespaces"""

    def __init__(self, ontology, namespace_domain, *, bel_function=None, encoding=None, ols_base=None, auth=None):
        """
        :param str ontology: The name of the ontology. Ex: ``uberon``, ``go``, etc.
        :param str namespace_domain: One of: :data:`pybel.constants.NAMESPACE_DOMAIN_BIOPROCESS`,
                            :data:`pybel.constants.NAMESPACE_DOMAIN_CHEMICAL`,
                            :data:`pybel.constants.NAMESPACE_DOMAIN_GENE`, or
                            :data:`pybel.constants.NAMESPACE_DOMAIN_OTHER`
        :param str bel_function: The BEL function of elements of this ontology. One of
                                :data:`pybel.constants.ABUNDANCE`, etc.
        :param str ols_base: An optional, custom OLS base url
        :param tuple[str,str] auth: A pair of (str username, str password) to give to the auth keyword of the
                                    constructor of :class:`artifactory.ArtifactoryPath`. Defaults to the result of
                                    :func:`pybel_tools.resources.get_arty_auth`.
        """
        if bel_function is None and encoding is None:
            raise ValueError('either bel_function or encoding must be specified')

        self.ontology = ontology
        self.namespace_domain = namespace_domain
        self._bel_function = bel_function
        self._encoding = encoding
        self.ols_client = OlsClient(ols_base=ols_base)
        self.auth = auth

        self.metadata = self.ols_client.get_ontology(self.ontology)

    @property
    def bel_function(self):
        return rev_abundance_labels[self._bel_function]

    @property
    def encodings(self):
        """The encodings that should be used when outputting as a BEL namespace"""
        if self._encoding:
            return self._encoding

        return function_to_encoding[self._bel_function]

    @property
    def title(self):
        """The ontology's full name"""
        return self.metadata['config']['title']

    @property
    def preferred_prefix(self):
        """The preferred prefix, usually uppercase"""
        return self.metadata['config']['preferredPrefix']

    @property
    def description(self):
        """A description of the ontology"""
        return self.metadata['config']['description']

    @property
    def version(self):
        """The version of the ontology"""
        return self.metadata['config']['version']

    @property
    def version_iri(self):
        """The IRI of this version of the ontoloy"""
        return self.metadata['config']['versionIri']

    def _get_values(self):
        return self.ols_client.iter_labels(self.ontology)

    def write_namespace(self, file=None):
        """Serializes a BEL namespace.

        :param file file: A write-enable file or file-like
        """
        values = self._get_values()

        write_namespace(
            namespace_name=self.title,
            namespace_keyword=self.preferred_prefix,
            namespace_domain=self.namespace_domain,
            namespace_description=self.description,
            namespace_version=self.version,
            citation_name=self.title,
            citation_url=self.version_iri,
            author_name='Charles Tapley Hoyt',
            author_contact='charles.hoyt@scai.fraunhofer.de',
            author_copyright='Creative Commons by 4.0',
            values=values,
            functions=self.encodings,
            file=file
        )

    def write_annotation(self, file=None):
        """Serializes a BEL annotation

        :param file file: A write-enable file or file-like
        """
        write_annotation(
            keyword=self.preferred_prefix,
            values={x: '' for x in self._get_values()},
            citation_name=self.title,
            description=self.description,
            version=self.version,
            author_name='Charles Tapley Hoyt',
            author_contact='charles.hoyt@scai.fraunhofer.de',
            author_copyright='Creative Commons by 4.0',
            file=file
        )

    def _write_hierarchy_header(self, file=None):
        write_boilerplate(
            name=self.title,
            description=self.description,
            authors='Charles Tapley Hoyt',
            contact='charles.hoyt@scai.fraunhofer.de',
            version=self.version,
            namespace_url={self.preferred_prefix: get_latest_arty_namespace(self.ontology)},
            annotation_url={},
            annotation_patterns={},
            namespace_patterns={},
            file=file
        )

    def _get_hierarchy(self):
        return self.ols_client.iter_hierarchy(self.ontology)

    def _write_hierarchy_body(self, file=None):
        print('SET Citation = {{"{}","{}"}}'.format(CITATION_TYPE_URL, self.version_iri), file=file)
        print('SET Evidence = "Automatically generated hierarchy from {}"\n'.format(self.version_iri), file=file)

        for parent, child in self._get_hierarchy():
            print(
                '{fn}({keyword}:{child}) {relation} {fn}({keyword}:{parent})'.format(
                    fn=self.bel_function,
                    keyword=self.preferred_prefix,
                    relation=IS_A,
                    parent=ensure_quotes(parent),
                    child=ensure_quotes(child)
                ),
                file=file
            )

    def write_hierarchy(self, file=None):
        """Serializes the hierarchy

        :param file file: A write-enable file or file-like
        """
        self._write_hierarchy_header(file=file)
        self._write_hierarchy_body(file=file)

    def deploy_namespace(self, hash_check=True):
        """Gets the data and writes BEL namespace file to Artifactory

        :param bool hash_check: Ensure the hash is unique before deploying
        :return: The path, if it was deployed successfully, else none.
        :rtype: str or None
        """
        file_name = get_today_arty_namespace(self.ontology)

        with open(file_name, 'w') as file:
            self.write_namespace(file)

        return deploy_namespace(file_name, self.ontology, hash_check=hash_check, auth=self.auth)

    def deploy_annotation(self, hash_check=True):
        """Gets the data and writes BEL annotation file to Artifactory

       :param bool hash_check: Ensure the hash is unique before deploying
       :return: The path, if it was deployed successfully, else none.
       :rtype: str or None
       """
        file_name = get_today_arty_annotation(self.ontology)

        with open(file_name, 'w') as file:
            self.write_annotation(file)

        return deploy_annotation(file_name, self.ontology, hash_check=hash_check, auth=self.auth)

    def deploy_hierarchy(self):
        """Gets the data and writes BEL hierarchy file to Artifactory

        :return: The path, if it was deployed successfully, else none.
        :rtype: str or None
        """
        file_name = get_today_arty_knowledge(self.ontology)

        with open(file_name, 'w') as file:
            self.write_hierarchy(file)

        return deploy_knowledge(file_name, self.ontology, auth=self.auth)

    def deploy(self, hash_check=True):
        """Gets the data and writes BEL namespace file to Artifactory

        :param bool hash_check: Ensure the hash is unique before deploying
        """
        namespace_url = self.deploy_namespace(hash_check=hash_check)

        if not namespace_url:
            self.deploy_hierarchy()


class OlsConstrainedNamespaceOntology(OlsNamespaceOntology):
    """Specifies that only the hierarchy under a certain term should be followed"""

    def __init__(self, ontology, namespace_domain, base_term_iri, *, bel_function=None, encoding=None, ols_base=None,
                 auth=None):
        """
        :param str ontology: The name of the ontology. Ex: ``uberon``, ``go``, etc.
        :param str namespace_domain: One of: :data:`pybel.constants.NAMESPACE_DOMAIN_BIOPROCESS`,
                            :data:`pybel.constants.NAMESPACE_DOMAIN_CHEMICAL`,
                            :data:`pybel.constants.NAMESPACE_DOMAIN_GENE`, or
                            :data:`pybel.constants.NAMESPACE_DOMAIN_OTHER`
        :param str bel_function: The BEL function of elements of this ontology
        :param str ols_base: An optional, custom OLS base url
        :param tuple[str,str] auth: A pair of (str username, str password) to give to the auth keyword of the
                                    constructor of :class:`artifactory.ArtifactoryPath`. Defaults to the result of
                                    :func:`pybel_tools.resources.get_arty_auth`.
        """
        super(OlsConstrainedNamespaceOntology, self).__init__(
            ontology=ontology,
            namespace_domain=namespace_domain,
            bel_function=bel_function,
            encoding=encoding,
            ols_base=ols_base,
            auth=auth
        )
        self.base_term_iri = base_term_iri

        self._labels = None
        self._hierarchy = None

    def _get_values(self):
        return self.ols_client.iter_descendants_labels(self.ontology, self.base_term_iri)

    def _get_hierarchy(self):
        raise NotImplementedError
