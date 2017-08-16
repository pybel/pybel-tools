# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import CITATION, CITATION_TYPE, CITATION_TYPE_PUBMED, CITATION_REFERENCE
from pybel.manager.models import Citation
from pybel_tools.citation_utils import get_citations_by_pmids
from tests.constants import ManagerMixin


class TestCitations(ManagerMixin):
    def test_enrich(self):
        g = BELGraph()

        g.add_node(1)
        g.add_node(2)

        g.add_edge(1, 2, attr_dict={
            CITATION: {
                CITATION_TYPE: CITATION_TYPE_PUBMED,
                CITATION_REFERENCE: "9611787"
            }
        })

        pmids = ["9611787"]

        stored_citations = self.manager.session.query(Citation).all()

        self.assertEqual(0, len(stored_citations))

        get_citations_by_pmids(pmids, manager=self.manager)

        stored_citations = self.manager.session.query(Citation).all()
        self.assertEqual(1, len(stored_citations))


if __name__ == '__main__':
    unittest.main()
