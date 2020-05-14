# -*- coding: utf-8 -*-

"""Tests for the HTML summary assembler."""

import os
import tempfile
import unittest

from pybel.examples import sialic_acid_graph

try:
    import bio2bel_hgnc
    import bio2bel_entrez
except ImportError:
    bio2bel_hgnc = None
    bio2bel_entrez = None


class TestSummaryAssembler(unittest.TestCase):
    """Tests for the assemblers."""

    def test_summary_to_html_path(self):
        """Test to_html_path."""
        import pybel_tools.assembler.html

        with tempfile.TemporaryDirectory() as tmpdirname:
            path = os.path.join(tmpdirname, 'summary.html')
            pybel_tools.assembler.html.to_html_path(graph=sialic_acid_graph, path=path)

            self.assertTrue(os.path.exists(path))
            with open(path) as file:
                contents = file.read()
                self.assertIn('<html', contents)
                self.assertIn('PTPN6', contents)

    @unittest.skip('Need to upgrade Bio2BEL to PyBEL 14')
    def test_ideogram_to_html_path(self):
        """Test to_html_path."""
        import pybel_tools.assembler.ideogram

        with tempfile.TemporaryDirectory() as tmpdirname:
            path = os.path.join(tmpdirname, 'ideogram.html')
            pybel_tools.assembler.ideogram.to_html_path(graph=sialic_acid_graph, path=path)

            self.assertTrue(os.path.exists(path))
            with open(path) as file:
                contents = file.read()
                self.assertIn('<html', contents)
                self.assertIn('PTPN6', contents)
