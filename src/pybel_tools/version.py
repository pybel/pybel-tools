
__all__ = [
    'VERSION',
    'get_version',
]

VERSION = '0.7.4-dev'

def get_version() -> str:
    """Get the current PyBEL Tools version."""
    return VERSION
