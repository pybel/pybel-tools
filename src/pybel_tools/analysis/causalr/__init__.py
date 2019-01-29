# -*- coding: utf-8 -*-

"""An implementation of the CausalR algorithm described by [Bradley2017]_.

.. [Bradley2017] Bradley, G., & Barrett, S. J. (2017). `CausalR - extracting mechanistic sense from genome scale
                 data <https://doi.org/10.1093/bioinformatics/btx425>`_. Bioinformatics, (June), 1–3.
"""

from .algorithm import rank_causalr_hypothesis

__all__ = [
    'rank_causalr_hypothesis',
]
