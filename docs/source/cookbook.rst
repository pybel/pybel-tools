Cookbook
========

Pickling lots of BEL scripts
----------------------------
All of the BEL scripts in the current working directory (and sub-directories) can be pickled in-place with the
following command (add :code:`-d` to specify a different directory)

.. code-block:: sh

    $ python3 -m pybel_tools io convert -d ~/bms/aetionomy/


Getting Data in to the Cache
----------------------------
Before running the service, some data can be pre-loaded in your cache.

Loading Selventa Corpra
~~~~~~~~~~~~~~~~~~~~~~~
The Selventa Small Corpus and Large Corpus are two example BEL documents distributed by the
`OpenBEL framework <https://wiki.openbel.org/display/home/Summary+of+Large+and+Small+BEL+Corpuses>`_. They are good
examples of many types of BEL statements and can be used immediately to begin exploring. Add :code:`-v` for more
logging information during compilation. This is highly suggested for the first run, since it takes a while to cache
all of the namespaces and annotations. This only has to be done once, and will be much faster the second time!

Small Corpus
************
.. code-block:: sh

    $ python3 -m pybel_tools ensure small_corpus -v

Large Corpus
************
.. code-block:: sh

    $ python3 -m pybel_tools ensure large_corpus -v

Loading Other Resources
~~~~~~~~~~~~~~~~~~~~~~~
Gene Families
*************
.. code-block:: sh

    $ python3 -m pybel_tools ensure gene_families -v

Named Protein Complexes
***********************
.. code-block:: sh

    $ python3 -m pybel_tools ensure named_complexes -v

Orthology Relations
*******************
Coming soon!

Uploading Precompiled BEL
~~~~~~~~~~~~~~~~~~~~~~~~~
A single network stored as a PyBEL gpickle can quickly be uploaded using the following code:

.. code-block:: sh

    $ python3 -m pybel_tools io upload -p /path/to/my_network.gpickle

Uploading Multiple Networks
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Multiple networks in a given directory and sub-directories can be uploaded by adding the :code:`-r` tag.

.. code::

    $ python3 -m pybel_tools io upload -p ~/bms/aetionomy/ -r