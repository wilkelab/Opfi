Getting Started
===============

.. _installation:

Installation
------------

The recommended way to install Opfi is with `Bioconda <https://bioconda.github.io/>`_, which requires the `conda <https://docs.conda.io/en/latest/>`_ package manager. This will install Opfi and all of its dependencies (which you can read more about below, see :ref:`dependencies`).

Currently, Bioconda supports only 64-bit Linux and Mac OS. Windows users can still install Opfi with pip (see below); however, the complete installation procedure has not been fully tested on a Windows system. 

.. _install-with-conda:

Install with conda (Linux and Mac OS only)
##########################################

First, set up conda and Bioconda following the `quickstart <https://bioconda.github.io/user/install.html>`_ guide. Once this is done, run:

.. code-block:: bash

    conda install -c bioconda opfi

And that's it! Note that this will install Opfi in the conda environment that is currently active. To create a fresh environment with Opfi installed, do:

.. code-block:: bash

    conda create --name opfi-env -c bioconda opfi
    conda activate opfi-env

.. _install-with-pip:

Install with pip
################

This method does not automatically install non-Python dependencies, so they will need to be installed separately, following their individual installation instructions. A complete list of required software is provided below, see :ref:`dependencies`. Once this step is complete, install Opfi with pip by running:

.. code-block:: bash

    pip install opfi

Install from source
###################

Finally, the latest development build may be installed directly from Github. First, non-Python :ref:`dependencies` will need to be installed in the working environment. An easy way to do this is to first install Opfi with conda using the :ref:`install-with-conda` method (we'll re-install the development version of the Opfi package in the next step). Alternatively, dependencies can be installed individually.

Once dependencies have been installed in the working environment, run the following code to download and install the development build:

.. code-block:: bash

    git clone https://github.com/wilkelab/Opfi.git
    cd Opfi
    pip install . # or pip install -e . for an editable version
    pip install -r requirements # if conda was used, this can be skipped

Testing the build
#################

Regardless of installation method, users can download and run Opfi's suite of unit tests to confirm that the build is working as expected. First download the tests from Github:

.. code-block:: bash

    git clone https://github.com/wilkelab/Opfi
    cd Opfi

And then run the test suite using pytest:

.. code-block:: bash

    pytest --runslow --runmmseqs --rundiamond

This may take a minute or so to complete. 

.. _dependencies:    

Dependencies
------------

Opfi uses the following bioinformatics software packages to find and annotate genomic features:

.. csv-table:: 
   :header: "Application", "Description"

   "`NCBI BLAST+ <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs>`_", "Protein and nucleic acid homology search tool"
   "`Diamond <https://github.com/bbuchfink/diamond>`_", "Alternative to BLAST+ for fast protein homology searches"
   "`MMseqs2 <https://github.com/soedinglab/MMseqs2>`_", "Alternative to BLAST+ for fast protein homology searches"
   "`PILER-CR <https://www.drive5.com/pilercr/>`_", "CRISPR repeat detection"
   "`Generic Repeat Finder <https://github.com/bioinfolabmu/GenericRepeatFinder>`_", "Transposon-associated repeat detection"

The first three (BLAST+, Diamond, and MMseqs2) are popular homology search tools, and are used in :class:`gene_finder.pipeline.Pipeline` for gene annotation. The user specifies which homology search tool to use during pipeline setup (see :class:`gene_finder.pipeline.Pipeline` for details). Note that the BLAST+ distribution contains multiple utilities for homology searches, three of which (blastp, blastn, and PSI-BLAST) are currently supported by Opfi. 

The primary difference between these tools is in their speed/sensitivity trade-offs. The following table is based on our own user experience, and may help users decide which tool best meets their needs. Note that performance tests are inherently hardware and context dependent, so this should be taken as a loose guide, rather than a definitive comparison. 

.. csv-table::
    :header: "Application", "Relative sensitivity", "Relative speed"

    "Diamond", "+", "++++"
    "MMseqs2", "++", "+++"
    "blastp", "+++", "++"
    "PSI-BLAST", "++++", "+"

The last two software dependencies, PILER-CR and Generic Repeat Finder (GRF), deal with annotation of repetive sequences in DNA. PILER-CR can identify CRISPR arrays, regions of alternatating ~30 bp direct repeat and variable sequences that play a role in prokaryotic immunity. GRF identifies repetitive sequences..
