# -*- coding: utf-8 -*-

"""An implementation of chemical-based mechanism enrichment with NeuroMMSig described by [Hoyt2018]_.

This algorithm has multiple steps:

1. Select NeuroMMSig networks for AD, PD, and epilepsy
2. Select drugs from DrugBank, and their targets
3. Run NeuroMMSig algorithm on target list for each network and each mechanism
4. Store in database

.. [Hoyt2018] Hoyt, C. T., *et al*. (2018) `A systematic approach for identifying shared mechanisms in epilepsy and its
              comorbidities <https://doi.org/10.1093/database/bay050>`_, Database, Volume 2018, 1 January 2018, bay050
"""

from .algorithm import multi_run_epicom, run_epicom

__all__ = [
    'run_epicom',
    'multi_run_epicom',
]
