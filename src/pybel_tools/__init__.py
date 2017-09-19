# -*- coding: utf-8 -*-

"""

PyBEL Tools is tested on Python3 installations on Mac OS and Linux on 
`Travis CI <https://travis-ci.org/pybel/pybel-tools>`_.

.. warning:: Python2 and Windows are not thoroughly tested

Installation
------------

Easiest
~~~~~~~
Download the latest stable code from `PyPI <https://pypi.python.org/pypi/pybel-tools>`_ with:

.. code-block:: sh

   $ python3 -m pip install pybel_tools

Get the Latest
~~~~~~~~~~~~~~~
Download the most recent code from `GitHub <https://github.com/pybel/pybel-tools>`_ with:

.. code-block:: sh

   $ python3 -m pip install git+https://github.com/pybel/pybel-tools.git@develop
   
   
For Developers
~~~~~~~~~~~~~~
Clone the repository from `GitHub <https://github.com/pybel/pybel-tools>`_ and install in editable mode with:

.. code-block:: sh

   $ git clone https://github.com/pybel/pybel-tools.git@develop
   $ cd pybel-tools
   $ python3 -m pip install -e .
   
Caveats
-------
PyBEL Tools contains many dependencies, including the scientific Python Stack (numpy, scipy, etc.). This makes 
installation difficult for Windows users, for whom Python cannot easily build C extensions. We recommend using an 
`Anaconda <https://www.continuum.io/downloads>`_ distribution of Python, which includes these precompiled.

Adding New Web Services
-----------------------
These documents are for developers who want to extend PyBEL Web.

Admin Only
~~~~~~~~~~
For making web pages that can only be seen by certain user types, see: 
`roles required <https://pythonhosted.org/Flask-Security/api.html#flask_security.decorators.roles_required>`_ 
documentation on Flask-Security.
"""

from . import analysis
from . import api
from . import citation_utils
from . import comparison
from . import definition_utils
from . import document_utils
from . import filters
from . import generation
from . import integration
from . import ioutils
from . import mutation
from . import orthology
from . import query
from . import recuration
from . import selection
from . import serialization
from . import summary
from . import utils
from . import visualization

__all__ = (
    analysis,
    api,
    citation_utils,
    comparison,
    definition_utils,
    document_utils,
    filters,
    generation,
    integration,
    ioutils,
    mutation,
    orthology,
    query,
    recuration,
    selection,
    serialization,
    summary,
    utils,
    visualization,
)

__version__ = '0.4.0'

__title__ = 'pybel_tools'
__description__ = 'Tools for using BEL documents in python'
__url__ = 'https://github.com/pybel/pybel-tools'

__author__ = 'Charles Tapley Hoyt'
__email__ = 'cthoyt@gmail.com'

__license__ = 'Apache License 2.0'
__copyright__ = 'Copyright (c) 2016 Charles Tapley Hoyt'
