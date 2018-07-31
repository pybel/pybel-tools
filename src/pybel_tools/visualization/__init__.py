# -*- coding: utf-8 -*-

"""A wrapper around PyBEL-Jupyter."""

import warnings

__all__ = [
    'to_jupyter',
    'to_html',
]


def to_jupyter(*args, **kwargs):
    warnings.warn('use pybel-jupyter package instead', DeprecationWarning)
    from pybel_jupyter import to_jupyter
    return to_jupyter(*args, **kwargs)


def to_html(*args, **kwargs):
    warnings.warn('use pybel-jupyter package instead', DeprecationWarning)
    from pybel_jupyter import to_html
    return to_html(*args, **kwargs)
