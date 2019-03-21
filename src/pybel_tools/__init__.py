# -*- coding: utf-8 -*-

"""A cool pool of tools for PyBEL.

Installation
------------
Easiest
~~~~~~~
Download the latest stable code from `PyPI <https://pypi.python.org/pypi/pybel-tools>`_ with:

.. code-block:: sh

   $ pip install pybel_tools

Get the Latest
~~~~~~~~~~~~~~~
Download the most recent code from `GitHub <https://github.com/pybel/pybel-tools>`_ with:

.. code-block:: sh

   $ pip install git+https://github.com/pybel/pybel-tools.git


For Developers
~~~~~~~~~~~~~~
Clone the repository from `GitHub <https://github.com/pybel/pybel-tools>`_ and install in editable mode with:

.. code-block:: sh

   $ git clone https://github.com/pybel/pybel-tools.git
   $ cd pybel-tools
   $ pip install -e .

Caveats
~~~~~~~
PyBEL Tools contains many dependencies, including the scientific Python Stack (numpy, scipy, etc.). This makes
installation difficult for Windows users, for whom Python cannot easily build C extensions. We recommend using an
`Anaconda <https://www.continuum.io/downloads>`_ distribution of Python, which includes these precompiled.

Testing
-------
PyBEL-Tools is tested on Python3 on Linux on
`Travis CI <https://travis-ci.org/pybel/pybel-tools>`_.
"""

from .utils import get_version
