---
title: 'Opfi: A Python package for identifying gene clusters in large genomics and metagenomics data sets'
tags:
  - Python
  - bioinformatics
  - metagenomics
  - gene cluster analysis
authors:
  - name: Alexis M. Hill^[co-first author, corresponding author]
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
date: 28 July 2021
bibliography: paper.bib

---

# Summary

Gene clusters perform a diverse set of functions, many of which are relevant to biotechnology, but tools for identification of gene clusters are lacking. Therefore, we developed Opfi: a modular pipeline for identification of arbitrary gene clusters in assembled genomic or metagenomic data. Opfi contains functions for de-novo annotation, de-deduplication, and visualization of putative gene clusters. It utilizes a customizable rule-based filtering approach for selection of candidate systems that adhere to user-defined criteria. Opfi is implemented in Python, and is available on the Python Package Index and the Bioconda channel.  

# Statement of need

Gene clusters are sets of co-localized, often contiguous genes that together perform specific functions. These systems have been applied successfully to meet a number of biotechnological needs, including biofuel production, organic compound synthesis, and gene editing [@Fischbach:2010].  Although there are many tools available for annotation of singular genes (or protein domains) in biological sequence data [@Camacho:2009; @Steinegger:2017; @Buchfink:2021], these programs do not identify whole gene clusters out of the box. In many cases, researchers must combine bioinformatics tools ad-hoc, resulting in one-off pipelines that can be difficult to reproduce. Several software packages have been developed to discovery of specific gene clusters [@Blin:2019; @Santos-Aberturas:2019; @vanHeel:2018], but these tools may not be sufficiently flexible to identify clusters of an arbitrary genomic composition. To address these gaps, we developed a modular pipeline that integrates multiple bioinformatics tools, providing a flexible, uniform computational framework for identification of arbitrary gene clusters. In a recent study, we used Opfi to uncover novel CRISPR-associated transposons (CASTs) in a large metagenomics dataset [@Rybarski:2021].

# Implementation

Opfi is implemented in Python, and uses several bioinformatics tools for gene cluster feature annotation [Camacho:2009; Steinegger:2017; Buchfink:2021; Edgar:2007; Shi:2019]. Users can install Opfi and all of it's dependencies from Bioconda, using the conda package manager. Opfi consists of two major components: Gene Finder, for discovery of gene clusters, and Operon Analyzer, for rule-based filtering, deduplication, visualization, and re-annotation of gene clusters identified by Gene Finder. All modules generate output in a comma-separated (CSV) format that is common to the entire package.

## Example Gene Finder usage

The following example script searches for putative CRISPR-Cas loci in the genome of *Rippkaea orientalis PCC 8802*. Information about the biological significance of this example, as well as data inputs and descriptions, can be found in the `tutorials` directory in the project GitHub repository. The example illustrates the use of the `Pipeline` class for setting up a gene cluster search. First, `add_seed_step` specifies a step to annotate *cas1* genes, using protein BLAST (BLASTP) [@Camacho:2009] and a database of representative Cas1 protein sequences. 10,000 bp regions directly up- and downstream of each putative *cas1* gene are selected for further analysis, and all other regions are discarded. Next, `add_filter_step` adds a step to annotate candidate regions for additonal *cas* genes. Candidates that do not have at least one additional *cas* gene are discarded from the master list of putative systems. Finally, `add_crispr_step` adds a step to search remaining candidates for CRISPR arrays, i.e. regions of alternatating ~30 bp direct repeat and variable sequences, using the PILER-CR repeat finding software [@Edgar:2007]. 

```python
from gene_finder.pipeline import Pipeline
import os

genomic_data = "GCF_000024045.1_ASM2404v1_genomic.fna.gz"
job_id = "r_orientalis"

p = Pipeline()
p.add_seed_step(db="cas1", name="cas1", e_val=0.001, blast_type="PROT")
p.add_filter_step(db="cas_all", name="cas", e_val=0.001, blast_type="PROT")
p.add_crispr_step()

p.run(job_id=job_id, data=genomic_data, span=10000, gzip=True)
```

Running this code creates the CSV file `r_orientalis_results.csv`, which contains information about each system identified; in this example, that is two canonical CRISPR-Cas systems, and one locus with weak homology to *cas* genes. Each line in the file represents a single putative feature in a candidate locus. Features from the same candidate are grouped together in the CSV. Detailed information about the output format can be found in the Opfi documentation: (<https://opfi.readthedocs.io/>).

## Example Operon Analyzer usage

In the previous example, passing systems must meet the relatively permissive criterion of having at least one *cas1* gene co-localized with one additional *cas* gene. This is sufficient to identify CRISPR-Cas loci, but may also capture regions that do not contain functional CRISPR-Cas systems, but rather consist of open reading frames (ORFs) with weak homology to *cas* genes. These improbable systems could be eliminated during the homology search by making the match acceptance threshold more restrictive (i.e by decreasing the e-value), however, this could result in the loss of interesting, highly diverged systems. Therefore, we implemented a module that enables post- homology search filtering of candidate systems, using flexible rules that can be combined to create sophisticated elimination functions. This allows the user to first perform a broad homology search with permissive parameters, and then apply rules to cull unlikely candidates without losing interesting and/or novel systems. Additionally, rules may be useful for selecting candidates with a specific genomic composition for downstream analysis. It should be noted that the use of the word "operon" throughout this library is an artifact from early development of Opfi. At this time, Opfi does not predict whether a candidate system represents a true operon, that is, a set of genes under control of the same promoter.

Rule-based filtering is illustrated with the following example. The sample script takes the output generated by the previous example and reconstructs each system as an `Operon` object. Next, the `RuleSet` class is used to assess each candidate; here, passing systems must contain two cascade genes (*cas5* and *cas7*) no more than 1000 bp apart, and at least one *cas3* (effector) gene. For a complete list of rules, see the Opfi [documentation](https://opfi.readthedocs.io/). 

```python
from operon_analyzer import analyze, rules

rs = rules.RuleSet()
rs.contains_group(["cas5", "cas7"], max_gap_distance_bp = 1000)
rs.require("cas3"))

with open("r_orientalis_results.csv", "r") as input_csv:
    with open("filtered_output.csv", "w") as output_csv:
        analyze.evaluate_rules_and_reserialize(input_csv, rs, output_csv)
```

After running this code, the file `filtered_output.csv` contains only high-confidence type-I CRISPR-Cas systems (re-serialized to CSV format) that passed all rules in the rule set. 

## Candidate visualization

Opfi integrates the `DNAFeaturesViewer` package [@Zulkower:2020] to create gene diagrams of candidate systems. Each input system is visualized as a single PNG image. The sample script below reads in output from the previous example, and generates two gene diagram images, one for each CRISPR-Cas system present in *Rippkaea orientalis*. One image is provided for reference in \autoref{fig:operon}. 

```python
from operon_analyzer import load, rules, visualize

feature_colors = { "cas1": "lightblue",
                    "cas2": "seagreen",
                    "cas3": "gold",
                    "cas4": "springgreen",
                    "cas5": "darkred",
                    "cas6": "thistle",
                    "cas7": "coral",
                    "cas8": "red",
                    "cas9": "palegreen",
                    "cas10": "blue",
                    "cas11": "tan",
                    "cas12": "orange",
                    "cas13": "saddlebrown",
                    "CRISPR array": "purple"
                    }

fs = rules.FilterSet().pick_overlapping_features_by_bit_score(0.9)
with open("filtered_output.csv", "r") as operon_data:
    operons = [operon for operon in load.load_operons(operon_data)]
    for operon in operons:
      fs.evaluate(operon)
    visualize.plot_operons(operons, output_directory=".", \
      plot_ignored=False, feature_colors=feature_colors)
```

The `FilterSet` class is used to resolve features with sequences that overlap by more than 90%. Specifically, only the overlapping feature with the highest bitscore value (a quantity that describes the overall quality of an alignment) is rendered when `pick_overlapping_features_by_bit_score` is applied. Note that is not a requirement for candidate visualization, but can improve gene diagram clarity.

![One of two type-I CRISPR-Cas systems present in the genome of *Rippkaea orientalis PCC 8802* (cyanobacteria). Note that the ORF beginning at position ~2500 has homology with both *cas1* and *cas4*. These alignments have identical bitscores (i.e. the goodness of alignments is quivalent, using this metric), so both annotations appear in the diagram, even though `pick_overlapping_features_by_bit_score` was applied.\label{fig:operon}](operon_diagram.png)

# Acknowledgements

The authors would like to thank the staff of the Texas Advanced Computing Center for providing computational resources, and members of the Finkelstein and Wilke labs for helpful discussions. This work was supported by an NIGMS grant R01GM124141 (to I.J.F.), the Welch Foundation grant F-1808 (to I.J.F.), NIGMS grant R01 GM088344 (to C.O.W.), and the College of Natural Sciences Catalyst Award for seed funding.

# References
