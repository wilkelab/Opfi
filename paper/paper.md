---
title: 'Opfi: A Python package for identification of gene cassettes in large genomics and metagenomics data sets'
tags:
  - Python
  - bioinformatics
  - metagenomics
  - gene cluster analysis
authors:
  - name: Alexis M. Hill
    affiliation: 1
  - name: James R. Rybarski^[co-first author]
    affiliation: 2
  - name: Kuang Hu
    affiliation: "1,2"
  - name: Ilya J. Finkelstein
    affiliation: "2,3"
  - name: Claus O. Wilke
    affiliation: 1
affiliations:
 - name: Department of Integrative Biology, The University of Texas at Austin, Austin, Texas 78712, USA
   index: 1
 - name: Department of Molecular Biosciences, The University of Texas at Austin, Austin, Texas 78712, USA
   index: 2
 - name: Center for Systems and Synthetic Biology, The University of Texas at Austin, Austin, Texas, 78712, USA
   index: 3
date: 15 June 2021
bibliography: paper.bib

---

# Summary

Advances in genome sequencing technology have led to an explosion of available metagenomics data, much of which may encode novel gene systems of academic and biotechnological interest. There is a critical need for tools that can systematically and reproducibly identify gene systems in large, often low-quality metagenomics data. Therefore, we developed Opfi: a modular, rule-based pipeline for identification of functional sets of genes, such as biosynthetic gene clusters or CRISPR-cas defense systems, in large genomics or metagenomics data sets. 

# Statement of need

Metagenomics is the sequencing and analysis of bacterial genomes sampled directly from the environment [@Thomas:2012]. This data often captures diversity not present in conventional lab-cultured genomics data, and is thus a rich source of novel genetic material with potentially useful properties [@Tringe:2005]. Gene cassettes, that is, sets of co-localized genes that together perform specific functions, are of particular interest to metagenomics surveyors. These systems can, for example, be leveraged to produce antimicrobial compounds or to facilitate gene editing, to name just a few applications. There are several homology search tools [@Camacho:2009; @Steinegger:2017; @Buchfink:2021] available for annotatation of singular genes (or protein domains), but these programs do not identify whole gene cassettes out of the box. In many cases, researchers must combine bioinformatics tools ad-hoc, resulting in one-off pipelines that can be difficult to reproduce. Several software packages have been developed to search for specific gene clusters [@Blin:2019; @Santos-Aberturas:2019; @vanHeel:2018], but these tools may not be sufficiently flexible to identify clusters of an arbitrary genomic composition. To address these gaps, we developed a modular pipeline that integrates multiple bioinformatics tools, providing a flexible, uniform computational framework for metagenomics discovery of arbitrary gene cassettes.

# Implementation

Opfi is implemented entirely in Python, and can be downloaded from the Python package index and installed directly in the user’s compute environment. It consists of two major components: Gene Finder, for discovery of novel gene cassettes, and Operon Analyzer, for rule-based filtering, deduplication, visualization, and re-annotation of gene cassettes identified by Gene Finder. All modules generate output in a comma-separated format that is common to the entire package.

## Example Gene Finder usage

The following example script searches for putative CRISPR-Cas loci in the Rippkaea orientalis PCC 8802 (cyanobacteria) genome. Information about the biological significance of this example, as well as data inputs and descriptions, can be found in the `tutorials` directory in the project GitHub repository. The example illustrates the use of the `Pipeline` class for setting up a gene cassette search. First, `add_seed_step` specifies a step to annotate cas1 genes, using NCBI Blastp and a database of representative Cas1 protein sequences. The regions directly up- and downstream of each Cas1 hit define the candidate search space, which is reasonable since Cas1 is a highly conserved protein common to the majority of known CRISPR-cas loci. Next, `add_filter_step` adds a step to annotate candidate regions for additonal cas genes. Regions that do not consist of at least one putative cas gene are discarded ("filtered") from the master list of putative systems. Finally, `add_crispr_step` adds a step to search remaining candidates for CRISPR arrays, i.e. regions of alternatating ~30 bp direct repeat and variable sequences, using the PILER-CR repeat finding software [@Edgar:2007]. 

```python
from gene_finder.pipeline import Pipeline
import os

genomic_data = "GCF_000024045.1_ASM2404v1_genomic.fna.gz"
job_id = "r_orientalis"

p = Pipeline()
p.add_seed_step(db="cas1", name="cas1", e_val=0.001, blast_type="PROT")
p.add_filter_step(db="cas_all", name="cas", e_val=0.001, blast_type="PROT")
p.add_crispr_step()

results = p.run(job_id=job_id, data=genomic_data, gzip=True)
```

## Example Operon Analyzer usage

In the previous example, passing systems must meet the relatively permissive criterion of having at least one cas1 gene co-localized with one additional cas gene. This is sufficient to identify CRISPR-Cas loci, but may capture regions that do not contain functional CRISPR-Cas systems. Improbale systems could be eliminated at the homology search phase by making the acceptance threshold more restrictive, or by reducing the size of the search space. However, we have implemented a module that facilitates rule-based filtering of candidate systems identified by Gene Finder, to reduce "false positive" candidates or select only those with a desired genomic oragnization. 

The following script takes the output generated by the previous example and reconstructs each systems as an `Operon` object. The `RuleSet` class is used to assess each system; here, passing candidates must contain two cascade genes (cas5 and cas7) no more than 1000 base pairs apart, and at least one cas3 (effector) gene. A complete list of rules can be found on the package documentation site (<https://opfi.readthedocs.io/>). In this example, only high-confidence type-I CRISPR-Cas systems are selected, and the passing systems are re-serialized to CSV format. 

```python
from operon_analyzer import analyze, rules

rs = rules.RuleSet()
rs.contains_group(["cas5", "cas7"], max_gap_distance_bp = 1000)
rs.require("cas3"))

with open("r_orientalis_results.csv", "r") as input_csv:
    with open("filtered_output.csv", "w") as output_csv:
        analyze.evaluate_rules_and_reserialize(input_csv, rs, output_csv)
```

Opfi integrates the `DNAFeaturesViewer` package [@Zulkower:2020] to generate gene diagrams of candidate systems. Each input system is visualized as a single PNG image. The sample script below reads in output from the previous example, and generates two gene diagram images, one for each CRISPR-Cas system present in Rippkaea orientalis. One image is provided for reference in \autoref{fig:operon}.

```python
from operon_analyzer import load, visualize

out_dir = "."
with open("filtered_output.csv", "r") as operon_data:
    operons = load.load_operons(operon_data)
    visualize.plot_operons(operons, out_dir, nucl_per_line=25000)
```

![One of two type-I CRISPR-Cas gene systems present in the genome of Rippkaea orientalis PCC 8802 (cyanobacteria).\label{fig:operon}](operon_diagram.png)

# Acknowledgements

The authors would like to thank the staff of the Texas Advanced Computing Center for providing computational resources, and members of the Finkelstein and Wilke labs for helpful discussions. This work was supported by an NIGMS grant R01GM124141 (to I.J.F.), the Welch Foundation grant F-1808 (to I.J.F.), NIGMS grant R01 GM088344 (to C.O.W.), and the College of Natural Sciences Catalyst Award for seed funding.

# References