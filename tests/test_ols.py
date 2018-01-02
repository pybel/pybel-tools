import unittest

from pybel.constants import NAMESPACE_DOMAIN_OTHER
from pybel_tools.ols_utils import OlsNamespaceOntology


class TestOls(unittest.TestCase):
    def test_import_ols(self):
        ns = OlsNamespaceOntology('uberon', NAMESPACE_DOMAIN_OTHER, encoding='A')
