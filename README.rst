PyBEL Tools |zenodo|
====================
`PyBEL Tools <https://pybel-tools.readthedocs.io/>`_ is a suite of tools built on top of
`PyBEL <https://pybel.readthedocs.io>`_ to facilitate data management, integration, and analysis. For examples,
see the `PyBEL Notebooks <https://github.com/pybel/pybel-notebooks>`_ repository.

=========== =============== ================== =======================
Stable      |stable_build|  |stable_coverage|  |stable_documentation|
Development |develop_build| |develop_coverage| |develop_documentation|
=========== =============== ================== =======================

Citation
--------
If you use PyBEL in your work, please cite:

.. [1] Hoyt, C. T., *et al.* (2017). `PyBEL: a Computational Framework for Biological Expression Language <https://doi.org/10.1093/bioinformatics/btx660>`_. Bioinformatics, 34(December), 1–2.

Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------
PyBEL Tools can be installed easily from PyPI_ with the following code in
your favorite terminal:

.. code-block:: sh

    $ python3 -m pip install pybel_tools

or from the latest code on `GitHub <https://github.com/pybel/pybel-tools>`_ with:

.. code-block:: sh

    $ python3 -m pip install git+https://github.com/pybel/pybel-tools.git@develop

See the `installation documentation <http://pybel-tools.readthedocs.io/en/stable/installation.html>`_ for more advanced
instructions.

Documentation and Examples
--------------------------
- Documentation at https://pybel-tools.readthedocs.io
- Cookbook at https://github.com/pybel/pybel-notebooks

Acknowledgements
----------------
This package was originally developed as part of the master's work of
`Charles Tapley Hoyt <https://github.com/cthoyt>`_ at `Fraunhofer SCAI <https://www.scai.fraunhofer.de/>`_.

Links
-----
- Documented on `Read the Docs <https://pybel-tools.readthedocs.io/>`_
- Versioned on `GitHub <https://github.com/pybel/pybel-tools>`_
- Tested on `Travis CI <https://travis-ci.org/pybel/pybel-tools>`_
- Distributed by PyPI_
- Chat on `Gitter <https://gitter.im/pybel/Lobby>`_


.. _PyPI: https://pypi.python.org/pypi/pybel_tools

.. |stable_build| image:: https://travis-ci.org/pybel/pybel-tools.svg?branch=master
    :target: https://travis-ci.org/pybel/pybel-tools
    :alt: Stable Build Status

.. |develop_build| image:: https://travis-ci.org/pybel/pybel-tools.svg?branch=develop
    :target: https://travis-ci.org/pybel/pybel-tools
    :alt: Development Build Status

.. |stable_coverage| image:: https://codecov.io/gh/pybel/pybel-tools/coverage.svg?branch=master
    :target: https://codecov.io/gh/pybel/pybel-tools?branch=master
    :alt: Stable Coverage Status

.. |develop_coverage| image:: https://codecov.io/gh/pybel/pybel-tools/coverage.svg?branch=develop
    :target: https://codecov.io/gh/pybel/pybel-tools?branch=develop
    :alt: Development Coverage Status

.. |stable_documentation| image:: https://readthedocs.org/projects/pybel-tools/badge/?version=stable
    :target: http://pybel-tools.readthedocs.io/en/stable/
    :alt: Stable Documentation Status

.. |develop_documentation| image:: https://readthedocs.org/projects/pybel-tools/badge/?version=latest
    :target: http://pybel-tools.readthedocs.io/en/latest/
    :alt: Development Documentation Status

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/pybel-tools.svg
    :alt: Stable Supported Python Versions

.. |pypi_version| image:: https://img.shields.io/pypi/v/pybel-tools.svg
    :alt: Current version on PyPI

.. |pypi_license| image:: https://img.shields.io/pypi/l/pybel-tools.svg
    :alt: Apache 2.0 License

.. |zenodo| image:: https://zenodo.org/badge/70473008.svg
    :target: https://zenodo.org/badge/latestdoi/70473008