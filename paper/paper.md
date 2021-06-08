# Summary:

Advances in genome sequencing technology have led to an explosion of available metagenomics data, encoding potentially novel genetic systems of academic and biotechnological interest. There is a critical need for tools that can systematically and reproducibly identify genetic systems in large, often low-quality metagenomics data. To this end, we developed Opfi: a modular, rule-based pipeline for identification of functional sets of genes, such as biosynthetic gene clusters or CRISPR-cas defense systems, in large metagenomics datasets. 

# Statement of need:

Metagenomics is the sequencing and analysis of bacterial genomes sampled directly from the environment (ZZZ Thomas, Gilbert, and Meyer 2012). This data often captures diversity not present in conventional lab-cultured genomics data, and is thus a rich source of novel genetic material with potentially useful properties. Gene clusters, that is, sets of co-localized, usually non-homologous genes that together perform specific functions, are of particular interest to metagenomics surveyors. These systems can, for example, be leveraged to produce antimicrobial compounds or to facilitate gene editing, to name just a few applications. Many popular homology based (ZZZ blast, diamond, mmseqs refs) and probabilistic (ZZZ hmmer ref) tools exist for annotation of singular genes (or protein domains), but cannot perform multi- gene cluster identification out of the box. In many cases, researchers must combine bioinformatics tools ad-hoc, resulting in one-off pipelines that can be difficult to reproduce. Several software packages have been developed to search for specific gene clusters (ZZZ examples), but these tools may not be sufficiently flexible to identify clusters of an arbitrary genomic composition. To address these gaps, we developed a modular pipeline that integrates multiple bioinformatics tools, providing a flexible, uniform computational framework for metagenomics discovery of arbitrary genetic systems.

# Implementation:

Opfi is implemented entirely in Python, and can be downloaded from the Python package index and installed directly in the user’s compute environment. It consists of two major modules: Gene Finder, for discovery of candidate gene systems, and Operon Analyzer, for filtering, deduplication, visualization, and re-annotation of systems identified by Gene Finder. Each module (or sub-module) generates output in a comma-separated format that is common to the entire package.

## Gene finder:

## Operon Analyzer:

# Acknowledgements:

# References