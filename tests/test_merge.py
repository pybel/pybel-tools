# -*- coding: utf-8 -*-

import os
import tempfile
import unittest

from pybel_tools.document_utils import write_boilerplate, merge


class TestBoilerplate(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        for item in os.listdir(self.dir):
            os.remove(os.path.join(self.dir, item))
        os.rmdir(self.dir)

    def test_1(self):
        b1 = os.path.join(self.dir, 'boilerplate1.bel')
        b2 = os.path.join(self.dir, 'boilerplate2.bel')
        b3 = os.path.join(self.dir, 'bp_merge.bel')

        pmids_1 = [26209472, 26940375]
        pmids_2 = [26839406, 26612754]

        with open(b1, 'w') as f:
            write_boilerplate(
                b1,
                'BP1',
                'cthoyt+1@gmail.com',
                'Boilerplate Test Document 1',
                'Charles Tapley Hoyt',
                pmids=pmids_1,
                file=f
            )

        with open(b2, 'w') as f:
            write_boilerplate(
                b2,
                'BP2',
                'cthoyt+2@gmail.com',
                'Boilerplate Test Document 2',
                'Charles Tapley Hoyt',
                pmids=pmids_2,
                file=f
            )

        merge(b3, [b1, b2])
