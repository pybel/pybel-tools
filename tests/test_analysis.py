# -*- coding: utf-8 -*-

"""Tests to make sure the analysis packages import."""

import logging
import unittest

logger = logging.getLogger(__name__)


class TestImports(unittest.TestCase):
    """Test imports."""

    def test_import_epicom(self):
        from pybel_tools.analysis import epicom
        logger.debug('loaded %s', epicom)

    def test_import_causalr(self):
        from pybel_tools.analysis import causalr
        logger.debug('loaded %s', causalr)

    def test_import_concordance(self):
        from pybel_tools.analysis import concordance
        logger.debug('loaded %s', concordance)

    def test_import_heat(self):
        from pybel_tools.analysis import heat
        logger.debug('loaded %s', heat)

    def test_import_mechanisms(self):
        from pybel_tools.analysis import mechanisms
        logger.debug('loaded %s', mechanisms)

    def test_import_rcr(self):
        from pybel_tools.analysis import rcr
        logger.debug('loaded %s', rcr)
