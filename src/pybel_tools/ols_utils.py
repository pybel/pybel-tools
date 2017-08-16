# -*- coding: utf-8 -*-

from ols_client import OlsClient

from .definition_utils import write_namespace
from .resources import get_today_arty_namespace, deploy_namespace


class OlsNamespaceOntology:
    """Wraps the functions needed to use the OLS to generate and deploy BEL namespaces"""

    def __init__(self, ontology_name, namespace_domain, functions=None, ols_base=None, auth=None):
        """
        :param str ontology_name: The name of the ontology. Ex: ``uberon``
        :param str namespace_domain: One of: :data:`pybel.constants.NAMESPACE_DOMAIN_BIOPROCESS`,
                            :data:`pybel.constants.NAMESPACE_DOMAIN_CHEMICAL`,
                            :data:`pybel.constants.NAMESPACE_DOMAIN_GENE`, or
                            :data:`pybel.constants.NAMESPACE_DOMAIN_OTHER`
        :param str functions: The encoding for the elements in this namespace. See
                              :data:`pybel.constants.belns_encodings`
        :param str ols_base: An optional, custom OLS base url
        :param tuple[str] auth: A pair of (str username, str password) to give to the auth keyword of the constructor of
                            :class:`artifactory.ArtifactoryPath`. Defaults to the result of
                            :func:`pybel_tools.resources.get_arty_auth`.
        """
        self.ontology_name = ontology_name
        self.namespace_domain = namespace_domain
        self.functions = functions
        self.ols_client = OlsClient(ols_base=ols_base)
        self.auth = auth

    def write(self, file):
        """Serializes a BEL namespace.

        :param file file: A write-enable file or file-like
        """
        metadata = self.ols_client.get_metadata(self.ontology_name)
        values = self.ols_client.get_labels(self.ontology_name)

        config = metadata['config']

        write_namespace(
            namespace_name=config['title'],
            namespace_keyword=config['preferredPrefix'],
            namespace_domain=self.namespace_domain,
            namespace_description=config['description'],
            namespace_version=config['version'],
            citation_name=config['title'],
            citation_url=config['versionIri'],
            author_name='Charles Tapley Hoyt',
            author_contact='charles.hoyt@scai.fraunhofer.de',
            author_copyright='Creative Commons by 4.0',
            values=values,
            functions=self.functions,
            file=file
        )

    def deploy(self, hash_check=True):
        """Gets the data and writes BEL namespace file to Artifactory

        :param bool hash_check: Ensure the hash is unique before deploying
        :return: The path, if it was deployed successfully, else none.
        :rtype: str
        """
        file_name = get_today_arty_namespace(self.ontology_name)

        with open(file_name, 'w') as file:
            self.write(file)

        return deploy_namespace(file_name, self.ontology_name, hash_check=hash_check, auth=self.auth)
