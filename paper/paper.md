---
title: 'Opfi: A Python package for constructing reproducible genomics and metagenomics mining pipelines'
tags:
  - Python
  - genome mining
  - metagenomics
  - gene cluster analysis
authors:
  - name: Alexis M. Hill^[co-first author]
    affiliation: 1
  - name: James R. Rybarski^[co-first author]
    affiliation: 2
  - name: Kuang Hu^[co-first author]
    affiliation: "1,2"
  - name: Ilya J. Finkelstein
    affiliation: "2,3"
  - name: Claus O. Wilke^[corresponding author]
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

# Summary:

Advances in genome sequencing technology have led to an explosion of available metagenomics data, encoding potentially novel genetic systems of academic and biotechnological interest. There is a critical need for tools that can systematically and reproducibly identify genetic systems in large, often low-quality metagenomics data. To this end, we developed Opfi: a modular, rule-based pipeline for identification of functional sets of genes, such as biosynthetic gene clusters or CRISPR-cas defense systems, in large metagenomics datasets. 

# Statement of need:

Metagenomics is the sequencing and analysis of bacterial genomes sampled directly from the environment (ZZZ Thomas, Gilbert, and Meyer 2012). This data often captures diversity not present in conventional lab-cultured genomics data, and is thus a rich source of novel genetic material with potentially useful properties. Gene clusters, that is, sets of co-localized, usually non-homologous genes that together perform specific functions, are of particular interest to metagenomics surveyors. These systems can, for example, be leveraged to produce antimicrobial compounds or to facilitate gene editing, to name just a few applications. Many popular homology based (ZZZ blast, diamond, mmseqs refs) and probabilistic (ZZZ hmmer ref) tools exist for annotation of singular genes (or protein domains), but cannot perform multi- gene cluster identification out of the box. In many cases, researchers must combine bioinformatics tools ad-hoc, resulting in one-off pipelines that can be difficult to reproduce. Several software packages have been developed to search for specific gene clusters (ZZZ examples), but these tools may not be sufficiently flexible to identify clusters of an arbitrary genomic composition. To address these gaps, we developed a modular pipeline that integrates multiple bioinformatics tools, providing a flexible, uniform computational framework for metagenomics discovery of arbitrary genetic systems.

# Implementation:

Opfi is implemented entirely in Python, and can be downloaded from the Python package index and installed directly in the user’s compute environment. It consists of two major modules: Gene Finder, for discovery of candidate gene systems, and Operon Analyzer, for rule-based filtering, deduplication, visualization, and re-annotation of systems identified by Gene Finder. Each module (or sub-module) generates output in a comma-separated format that is common to the entire package.

## Example Gene Finder usage:

```python
from gene_finder.pipeline import Pipeline

p = Pipeline()
p.add_seed_step(db="databases/common_genes", name="common", e_val=0.001, blast_type="PROT")
p.add_filter_step(db="databases/signature_genes", name="required", e_val=0.001, blast_type="PROT", min_prot_count=3)
p.add_blast_step(db="databases/ancillary_genes", name="ancillary", e_val=0.001, blast_type="PROT")

results = p.run(data="my_genome.fasta", output_directory="output")
```

## Example Operon Analyzer usage:

A toy example of rule-based filtering using Operon Analyzer:
```python
from operon_analyzer import analyze, rules

# Selects systems that adhere to the following conditions:
# 1. Must have an operon composed of gene1, gene2, and signature_gene_1
# 2. signature_gene_1 must be at least 3000 nucleotides in length
# 3. The system must also contain signature_gene_2, but it does not have to be in any particular position
# relative to the other required genes
rs = RuleSet().contains_group(feature_names = ["gene1", "gene2, signature_gene_1"], max_gap_distance_bp = 50) \
              .minimum_size("signature_gene_1", 3000))
              .require("signature_gene_2")

with open("gene_finder_output.csv", "r") as input_csv:
    with open(f"filtered_gene_finder_output.csv", "w") as output_csv:
        analyze.evaluate_rules_and_reserialize(input_csv, rs, output_csv)
```

Gene diagrams can be created for easy visualization of interesting systems:
```python
from operon_analyzer import load, visualize

with open("filtered_gene_finder_output.csv", "r") as operon_data:
    operons = load.load_operons(operon_data)
    visualize.plot_operons(operons=operons, output_directory="output", nucl_per_line=25000)
```

# Acknowledgements:

The authors would like to thank the staff of the Texas Advanced Computing Center for providing
computational resources, and members of the Finkelstein and Wilke labs for helpful
discussions. This work was supported by an NIGMS grant R01GM124141 (to I.J.F.), the Welch
Foundation grant F-1808 (to I.J.F.), NIGMS grant R01 GM088344 (to C.O.W.), and the College
of Natural Sciences Catalyst Award for seed funding.

# References