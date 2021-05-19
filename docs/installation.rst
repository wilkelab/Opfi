Getting Started
===============

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
    pip3 install -r requirements.txt

Dependencies
------------

Opfi makes use of several third-party softwares for finding and annotating genomic features. Depending on your use case, you may not need to install all of these; however, at a minimum users should have the NCBI BLAST+ application installed in their environment. The following table provides more details about required/optional dependencies, including links to associated application homepages.

.. csv-table:: 
   :header: "Application", "Required", "Description", "Anaconda distribution"

   "`NCBI BLAST+ <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs>`_", "Yes", "Protein and nucleic acid homology search tool", https://anaconda.org/bioconda/blast
   "`Diamond <https://github.com/bbuchfink/diamond>`_", "No", "Alternative to BLAST+ for fast protein homology searches", https://anaconda.org/bioconda/diamond
   "`MMseqs2 <https://github.com/soedinglab/MMseqs2>`_", "No", "Alternative to BLAST+ for fast protein homology searches", https://anaconda.org/bioconda/mmseqs2
   "`PILER-CR <https://www.drive5.com/pilercr/>`_", "No", "CRISPR repeat detection", https://anaconda.org/bioconda/piler-cr
   "`GenericRepeatFinder <https://github.com/bioinfolabmu/GenericRepeatFinder>`_", "No", "Transposon-associated repeat detection", "NA"

Testing your build
------------------

Users who opt to build Opfi from source can test their build by running the command `pytest --runslow` from the project root directory. Note that this assumes NCBI BLAST+ has already been installed. The following flags can direct pytest to test any optional depencies, as needed.

* `--rundiamond`: run tests that require Diamond.
* `--runmmseqs`: run tests that require MMseqs2.
* `--runpiler`: run tests that require PILER-CR.
* `--rungrf`: run tests that require GenericRepeatFinder.
