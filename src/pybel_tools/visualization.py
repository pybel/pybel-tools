# -*- coding: utf-8 -*-

"""A wrapper around PyBEL-Jupyter.

This functionality has been moved to: https://github.com/pybel/pybel-jupyter.
"""

import warnings

__all__ = [
    'to_jupyter',
    'to_html',
]


def to_jupyter(*args, **kwargs):
    warnings.warn('use pybel-jupyter package instead', DeprecationWarning)
    try:
        from pybel_jupyter import to_jupyter
    except ImportError:
        raise ImportError("""pybel_tools.visualization.to_jupyter has been moved to the pybel_jupyter module.

        Please install pybel_jupyter with pip install pybel-jupyter.
        """)
    else:
        return to_jupyter(*args, **kwargs)


def to_html(*args, **kwargs):
    warnings.warn('use pybel-jupyter package instead', DeprecationWarning)
    try:
        from pybel_jupyter import to_html
    except ImportError:
        raise ImportError("""pybel_tools.visualization.to_jupyter has been moved to the pybel_jupyter module.

        Please install pybel_jupyter with pip install pybel-jupyter.
        """)
    else:
        return to_html(*args, **kwargs)
