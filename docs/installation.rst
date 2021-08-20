Getting Started
===============

.. _installation:

Installation
------------

You can install Opfi with Pip:

.. code-block:: bash

    pip3 install opfi

Alternatively, you can install the latest version on Github:

.. code-block:: bash

    git clone https://github.com/alexismhill3/Opfi.git
    cd Opfi
    pip3 install .

Dependencies
------------

Opfi makes use of several third-party softwares for finding and annotating genomic features. Depending on your use case, you may not need to install all of these; however, at a minimum users should have the NCBI BLAST+ application installed in their environment. The following table provides more details about required/optional dependencies, including links to application homepages.

.. csv-table:: 
   :header: "Application", "Required", "Description", "Anaconda distribution"

   "`NCBI BLAST+ <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs>`_", "Yes", "Protein and nucleic acid homology search tool", https://anaconda.org/bioconda/blast
   "`Diamond <https://github.com/bbuchfink/diamond>`_", "No", "Alternative to BLAST+ for fast protein homology searches", https://anaconda.org/bioconda/diamond
   "`MMseqs2 <https://github.com/soedinglab/MMseqs2>`_", "No", "Alternative to BLAST+ for fast protein homology searches", https://anaconda.org/bioconda/mmseqs2
   "`PILER-CR <https://www.drive5.com/pilercr/>`_", "No", "CRISPR repeat detection", https://anaconda.org/bioconda/piler-cr
   "`GenericRepeatFinder <https://github.com/bioinfolabmu/GenericRepeatFinder>`_", "No", "Transposon-associated repeat detection", "NA"

Testing your build
------------------

Users who opt to build Opfi from source can test their build by running ``pytest`` from the project root directory. The following flags will direct pytest to run specific sets of tests (in addition to the core suite):

* ``--runmmseqs``: Run tests that require MMseqs2.
* ``--rundiamond``: Run tests that require Diamond.
* ``--runslow``: Run integration/end-to-end tests.
* ``--runprop``: Run very slow property tests.

For most users, running ``pytest --runslow`` is recommended. 

.. note::

    Several tests in the core suite require :program:`BLAST`, :program:`PILER-CR`, and/or :program:`GenericRepeatFinder`. Running ``pytest`` without first installing these dependencies will cause these tests to fail. 
