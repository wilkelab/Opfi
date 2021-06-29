---
title: 'Opfi: A Python package for constructing reproducible genomics and metagenomics mining pipelines'
tags:
  - Python
  - genome mining
  - metagenomics
  - gene cluster analysis
authors:
  - name: Alexis M. Hill^[co-first author]^[corresponding author]
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

Advances in genome sequencing technology have led to an explosion of available metagenomics data, encoding potentially novel genetic systems of academic and biotechnological interest. There is a critical need for tools that can systematically and reproducibly identify genetic systems in large, often low-quality metagenomics data. To this end, we developed Opfi: a modular, rule-based pipeline for identification of functional sets of genes, such as biosynthetic gene clusters or CRISPR-cas defense systems, in large metagenomics datasets. 

# Statement of need

Metagenomics is the sequencing and analysis of bacterial genomes sampled directly from the environment (ZZZ Thomas, Gilbert, and Meyer 2012). This data often captures diversity not present in conventional lab-cultured genomics data, and is thus a rich source of novel genetic material with potentially useful properties. Gene clusters, that is, sets of co-localized, usually non-homologous genes that together perform specific functions, are of particular interest to metagenomics surveyors. These systems can, for example, be leveraged to produce antimicrobial compounds or to facilitate gene editing, to name just a few applications. Many popular homology based (ZZZ blast, diamond, mmseqs refs) and probabilistic (ZZZ hmmer ref) tools exist for annotation of singular genes (or protein domains), but cannot perform multi- gene cluster identification out of the box. In many cases, researchers must combine bioinformatics tools ad-hoc, resulting in one-off pipelines that can be difficult to reproduce. Several software packages have been developed to search for specific gene clusters (ZZZ examples), but these tools may not be sufficiently flexible to identify clusters of an arbitrary genomic composition. To address these gaps, we developed a modular pipeline that integrates multiple bioinformatics tools, providing a flexible, uniform computational framework for metagenomics discovery of arbitrary genetic systems.

# Implementation

Opfi is implemented entirely in Python, and can be downloaded from the Python package index and installed directly in the user’s compute environment. It consists of two major modules: Gene Finder, for discovery of novel gene cassettes, and Operon Analyzer, for rule-based filtering, deduplication, visualization, and re-annotation of systems identified by Gene Finder. Each module (or sub-module) generates output in a comma-separated format that is common to the entire package.

## Example Gene Finder usage

The following example script searches for putative CRISPR-Cas loci in the Rippkaea orientalis PCC 8802 (cyanobacteria) genome. Information about the biological significance of this example, as well as data inputs and descriptions, can be found at https://github.com/wilkelab/Opfi/tutorials. The example illustrates the use of the `Pipeline` class for setting up a gene cassette search. First, `add_seed_step` specifies a step to annotate cas1 genes, using NCBI BlastP and a database of representative Cas1 protein sequences. The regions directly up- and downstream of each Cas1 hit define the candidate search space, which is reasonable since Cas1 is a highly conserved protein common to the majority of known CRISPR-cas loci. Next, `add_filter_step` adds a step to annotate candidate regions for additonal cas genes. Regions that do not consist of at least one putative cas gene are discarded ("filtered") from the master list of putative systems. Finally, `add_crispr_step` adds a step to search remaining candidates for CRISPR arrays, i.e. regions of alternatating ~30 bp direct repeat and variable sequences, using the PILER-CR repeat finding software (ZZZ PILER ref). 

```python
from gene_finder.pipeline import Pipeline
import os

genomic_data = "GCF_000024045.1_ASM2404v1_genomic.fna.gz"

p = Pipeline()
p.add_seed_step(db="cas1", name="cas1", e_val=0.001, blast_type="PROT", num_threads=1)
p.add_filter_step(db="cas_all", name="cas_all", e_val=0.001, blast_type="PROT", num_threads=1)
p.add_crispr_step()

# use the input filename as the job id
job_id = os.path.basename(genomic_data)
results = p.run(job_id=job_id, data=genomic_data, min_prot_len=90, span=10000, gzip=True)
```

## Example Operon Analyzer usage

In the previous example, passing systems must meet the relatively permissive criterion of having at least one cas1 gene co-localized with one additional cas gene. This is sufficient to identify CRISPR-Cas loci, but may capture regions that do not contain functional CRISPR-Cas systems. Improbale systems could be eliminated at the homology search phase by requiring the presence of additional genes, or by reducing the size of the search space. However, we have implemented a module that facilitates rule-based filtering of candidate systems identified by Gene Finder, to reduce "false positive" candidates or select only those with a desired genomic oragnization. 

The following script takes the output generated by the previous example and reconstructs each systems as an `Operon` object. The `RuleSet` class is used to assess each system; here, passing candidates must contain two cascade genes (cas5 and cas7) no more than 1000 bp apart, and at least one cas3 (effector) gene. A complete list of rules can be found at (ZZZ documetation site URL). In this example, only high-confidence type I CRISPR-cas systems are selected, and the passing systems are re-serialized to CSV format. 

```python
from operon_analyzer import analyze, rules

rs = rules.RuleSet().contains_group(feature_names = ["cas5", "cas7"], max_gap_distance_bp = 1000) \
                    .require("cas3"))

with open("GCF_000024045.1_ASM2404v1_genomic.fna.gz_results.csv", "r") as input_csv:
    with open("filtered_gene_finder_output.csv", "w") as output_csv:
        analyze.evaluate_rules_and_reserialize(input_csv, rs, output_csv)
```

Opfi integrates the `DNAFeaturesViewer` tool (https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) to generate gene diagrams of candidate systems. Each input system is visualized as a single png image. Below is a simple script to visualize passing systems from the previous code snippet. A sample gene diagram is also provided. 

```python
from operon_analyzer import load, visualize

with open("filtered_gene_finder_output.csv", "r") as operon_data:
    operons = load.load_operons(operon_data)
    visualize.plot_operons(operons=operons, output_directory=".", nucl_per_line=25000)
```

![Figure 1](operon_diagram.png)

# Acknowledgements

The authors would like to thank the staff of the Texas Advanced Computing Center for providing
computational resources, and members of the Finkelstein and Wilke labs for helpful
discussions. This work was supported by an NIGMS grant R01GM124141 (to I.J.F.), the Welch
Foundation grant F-1808 (to I.J.F.), NIGMS grant R01 GM088344 (to C.O.W.), and the College
of Natural Sciences Catalyst Award for seed funding.

# References