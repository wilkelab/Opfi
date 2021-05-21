.. Opfi documentation master file, created by
    sphinx-quickstart on Mon Aug 10 12:29:45 2020.
    You can adapt this file completely to your liking, but it should at least
    contain the root `toctree` directive.

Opfi - [insert interesting tagline here]
========================================

Opfi (short for "operon finder") is a bioinformatics tool for identifying and analyzing genomic neighborhoods of interest in huge metagenomic datasets. Opfi consists of two modules that are intended to be used together to create a complete bioinformatics discovery pipeline. `Gene Finder` provides an easy to use Python interface for executing blast searches, and is used to extract potential regions of interest from raw metegenomic data. `Operon Analyzer` provides a suite of methods for filtering, clustering, and visualizing candidates identified by `Gene Finder`.

Contents
--------

.. toctree::
    :maxdepth: 1
    
    installation
    examples
    modules