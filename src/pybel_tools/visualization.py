# -*- coding: utf-8 -*-

"""A wrapper around PyBEL-Jupyter.

This functionality has been moved to: https://github.com/pybel/pybel-jupyter.
"""

import warnings

from pybel_jupyter import to_html, to_jupyter

__all__ = [
    'to_jupyter',
    'to_html',
]

warnings.warn('The pybel_tools.visualization module has been externalized to pybel_jupyter.', DeprecationWarning)
