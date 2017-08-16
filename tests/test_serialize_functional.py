# -*- coding: utf-8 -*-

"""This module tests the serialization for BELIEF"""

import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.parser import BelParser
from pybel_tools.serialization.functional import convert_for_belief


class TestFunctionalize(unittest.TestCase):
    """These examples are taken from pybel's test functions in tests/test_parse_bel.py"""

    @classmethod
    def setUpClass(cls):
        cls.graph = BELGraph()
        cls.parser = BelParser(cls.graph)

    def setUp(self):
        self.parser.clear()
        self.parser.control_parser.citation = {
            CITATION_TYPE: 'Other'
        }
        self.parser.control_parser.evidence = 'test evidence'

    def test_tokens(self):
        """Checks token handling"""
        self.maxDiff = None

        tokens = {
            SUBJECT: {
                MODIFIER: ACTIVITY,
                TARGET: {
                    FUNCTION: PROTEIN,
                    IDENTIFIER: {NAMESPACE: 'HGNC', NAME: 'HMGCR'}
                },
                EFFECT: {
                    NAME: 'cat',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
            },
            RELATION: 'rateLimitingStepOf',
            OBJECT: {
                FUNCTION: BIOPROCESS,
                IDENTIFIER: {NAMESPACE: 'GOBP', NAME: 'cholesterol biosynthetic process'}
            }
        }

        expected = {
            'subject': {
                'term': {
                    'fx': 'activity',
                    'arguments': [
                        {
                            'term': {
                                'fx': 'proteinAbundance',
                                'arguments': [
                                    {
                                        'parameter': {
                                            'ns': 'HGNC',
                                            'value': 'HMGCR'
                                        }
                                    }
                                ]
                            }
                        },
                        {
                            'term': {
                                'fx': 'molecularActivity',
                                'arguments': [
                                    {
                                        'parameter': {
                                            'ns': 'bel',
                                            'value': 'catalyticActivity'
                                        }
                                    }
                                ]
                            }
                        }
                    ]
                }
            },
            'relationship': 'rateLimitingStepOf',
            'object': {
                'term': {
                    'fx': 'biologicalProcess',
                    'arguments': [
                        {
                            'parameter': {
                                'ns': 'GOBP',
                                'value': 'cholesterol biosynthetic process'
                            }
                        }
                    ]
                }
            }
        }

        self.assertEqual(expected, convert_for_belief(tokens))

    def test_full(self):
        """Checks full output from BEL statement"""
        self.maxDiff = None

        line = 'a(CHEBI:"phosphatidyl-L-serine") increases pep(p(HGNC:F3))'

        output = {
            "subject": {
                "term": {
                    "fx": "abundance",
                    "arguments": [
                        {
                            "parameter": {
                                "ns": "CHEBI",
                                "value": "phosphatidyl-L-serine"
                            }
                        }
                    ]
                }
            },
            "relationship": "increases",
            "object": {
                "term": {
                    "fx": "activity",
                    "arguments": [
                        {
                            "term": {
                                "fx": "proteinAbundance",
                                "arguments": [
                                    {
                                        "parameter": {
                                            "ns": "HGNC",
                                            "value": "F3"
                                        }
                                    }
                                ]
                            }
                        },
                        {
                            'term': {
                                'fx': 'molecularActivity',
                                'arguments': [
                                    {
                                        'parameter': {
                                            'ns': 'bel',
                                            'value': 'peptidaseActivity'
                                        }
                                    }
                                ]
                            }
                        }
                    ]
                }
            }
        }
        tokens = self.parser.parseString(line)
        self.assertEqual(output, convert_for_belief(tokens))
