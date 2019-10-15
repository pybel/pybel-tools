# -*- coding: utf-8 -*-

"""Version information for PyBEL-Tools."""

__all__ = [
    'VERSION',
    'get_version',
]

VERSION = '0.7.5'


def get_version() -> str:
    """Get the current PyBEL Tools version."""
    return VERSION
