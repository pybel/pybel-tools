# -*- coding: utf-8 -*-

import unittest

from pybel_tools.utils import min_tanimoto_set_similarity


class TestMinSimilarity(unittest.TestCase):
    def test_empty(self):
        a = {1, 2}
        b = set()

        self.assertEqual(0.0, min_tanimoto_set_similarity(a, b))

    def test_double_empty(self):
        a = set()
        b = set()
        self.assertEqual(0.0, min_tanimoto_set_similarity(a, b))

    def test_partial_equal(self):
        a = {1, 2, 3}
        b = {1, 2}
        self.assertEqual(1.0, min_tanimoto_set_similarity(a, b))

    def test_full_equal(self):
        a = {1, 2}
        b = {1, 2}
        self.assertEqual(1.0, min_tanimoto_set_similarity(a, b))
