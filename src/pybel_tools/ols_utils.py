# -*- coding: utf-8 -*-

from collections import defaultdict

from ols_client import OlsClient
from pybel.constants import CITATION_TYPE_URL, belns_encodings, IS_A
from pybel.parser.language import rev_abundance_labels
from pybel.utils import ensure_quotes
from .definition_utils import write_namespace
from .document_utils import write_boilerplate
from .resources import (
    get_today_arty_namespace,
    deploy_namespace,
    get_latest_arty_namespace,
    get_today_arty_knowledge,
    deploy_knowledge,
)

function_to_encoding = defaultdict(list)
for encoding, functions in belns_encodings.items():
    for fn in functions:
        function_to_encoding[fn].append(encoding)


class OlsNamespaceOntology:
    """Wraps the functions needed to use the OLS to generate and deploy BEL namespaces"""

    def __init__(self, ontology, namespace_domain, bel_function, ols_base=None, auth=None):
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
        self.ontology = ontology
        self.namespace_domain = namespace_domain
        self.function = rev_abundance_labels[bel_function]
        self.encodings = function_to_encoding[bel_function]

        self.ols_client = OlsClient(ols_base=ols_base)
        self.auth = auth

        self.metadata = self.ols_client.get_ontology(self.ontology)

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

    def write_namespace(self, file=None):
        """Serializes a BEL namespace.

        :param file file: A write-enable file or file-like
        """
        values = self.ols_client.iter_labels(self.ontology)

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

    def _write_hierarchy_header(self, file=None):
        write_boilerplate(
            document_name=self.title,
            description=self.description,
            authors='Charles Tapley Hoyt',
            contact='charles.hoyt@scai.fraunhofer.de',
            version=self.version,
            namespace_dict={self.preferred_prefix: get_latest_arty_namespace(self.ontology)},
            annotations_dict={},
            annotations_patterns={},
            namespace_patterns={},
            file=file
        )

    def _write_hierarchy_body(self, file=None):
        print('SET Citation = {{"{}","{}"}}'.format(CITATION_TYPE_URL, self.version_iri), file=file)
        print('SET Evidence = "Automatically generated hierarchy from {}"\n'.format(self.version_iri), file=file)

        for parent, child in self.ols_client.iter_hierarchy(self.ontology):
            print(
                '{fn}({keyword}:{child}) {relation} {fn}({keyword}:{parent})'.format(
                    fn=self.function,
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
