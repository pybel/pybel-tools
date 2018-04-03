# -*- coding: utf-8 -*-

import logging
import unittest

from pybel.examples import sialic_acid_graph, statin_graph
from pybel_tools.selection.paths import get_random_path

log = logging.getLogger(__name__)


class TestRandomPath(unittest.TestCase):

    def test_get_random_path(self):
        """Just make sure it doesn't crash"""
        for graph in (sialic_acid_graph, statin_graph):
            for _ in range(100):
                get_random_path(graph)


if __name__ == '__main__':
    unittest.main()
