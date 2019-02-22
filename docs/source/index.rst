PyBEL-Tools Documentation
=========================
`PyBEL-Tools <http://pybel-tools.readthedocs.io/>`_ is a suite of tools built on top of
`PyBEL <http://pybel.readthedocs.io>`_ to facilitate data management, integration, and analysis. For further examples,
see the `PyBEL-Notebooks <https://github.com/pybel/pybel-notebooks>`_ repository.

Citation
--------
If you use PyBEL and PyBEL Tools in your work, please cite [Hoyt2017]_:

.. [Hoyt2017] Hoyt, C. T., *et al.* (2017). `PyBEL: a Computational Framework for Biological Expression Language
              <https://doi.org/10.1093/bioinformatics/btx660>`_. Bioinformatics, 34(December), 1–2.

Links
-----
- Documented on `Read the Docs <http://pybel-tools.readthedocs.io/>`_
- Versioned on `GitHub <https://github.com/pybel/pybel-tools>`_
- Tested on `Travis CI <https://travis-ci.org/pybel/pybel-tools>`_
- Distributed by PyPI_
- Chat on `Gitter <https://gitter.im/pybel/Lobby>`_

.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   :name: springboard

   installation
   cli

.. toctree::
   :maxdepth: 2
   :caption: Graph Tools
   :name: graphutils

   summary
   filters
   selection
   integration
   mutation

.. toctree::
   :caption: Algorithms
   :name: algorithms

   algorithms/rcr
   algorithms/causalr
   algorithms/spia
   algorithms/neurommsig
   algorithms/epicom

.. toctree::
   :caption: Workflows
   :name: workflows

   stability
   expansion
   generation
   heat

.. toctree::
   :caption: Assemblers
   :name: assemblers

   assemblers/html
   assemblers/ideogram

.. toctree::
   :caption: Other Utilities
   :name: utilities

   documentutils
   utilities


Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _PyPI: https://pypi.python.org/pypi/pybel_tools
