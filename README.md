# Opfi

## Requirements

At a minimum, the NCBI BLAST+ software suite should be installed and on the user's PATH. BLAST+ installation instruction can be found [here](https://www.ncbi.nlm.nih.gov/books/NBK279671/). For annotation of CRISPR arrays, Opfi uses PILER-CR, which can be downloaded from the software [home page](https://www.drive5.com/pilercr/). A modified version of PILER-CR that detects mini (two repeat) CRISPR arrays is also available, and can be built with GNU make after cloning or downloading Opfi:

```
cd lib/pilercr1.06
sudo make install
```

## Installation

To install Opfi, simply clone the repository and run pip from the project root directory:

```
git clone https://github.com/wilkelab/Opfi.git
cd Opfi
pip3 install -r requirements.txt
pip3 install .
```

## Gene Finder usage

The gene finder API is centered around the Pipeline class. Typical usage looks something like this:
```python
from gene_finder.pipeline import Pipeline

p = Pipeline()
p.add_seed_step(db="data/blast_databases/tnsAB", name="tnsAB", e_val=0.001, type="PROT")
p.add_filter_step(db="data/blast_databases/cas", name="cas", e_val=0.001, type="PROT", min_prot_count=2)
p.add_blast_step(db="data/blast_databases/tnsCD", name="tnsCD", e_val=0.001, type="PROT")
p.add_crispr_step()

results = p.run(data="my_genome.fasta", min_prot_len=30, span=10000, outfrmt="CSV", outfile="mygenome.csv")
```
The docstrings in pipeline.py provide more details about pipeline usage, parameters, etc.
Example scripts to run jobs can be found under extras/.

# Operon Analyzer

Operon analyzer filters pipeline results and identifies promising candidate operons according to a given set of criteria, and also contains some tools for visualizing candidates and performing basic statistics.

## Analysis

The analysis module provides tools to identify operons that conform to certain rules, such as requiring that they contain a certain gene, or that two genes are within a given distance of each other (the full list is given below). CSV output is written to stdout, which identifies the outcome of the analysis for each putative operon.

Rules defined with the `RuleSet` determine whether an operon should be considered a candidate for further analysis. 
Filters defined with the `FilterSet` help define which features to consider when evaluating rules. You might, for example, want to exclude any operon containing a particular gene, but if a copy of that gene coincidentally exists 5 kb from the true operon, you might want to ignore it for the purposes of evaluating your rules. 

A sample script that performs this task is given here:

```python
import sys
from operon_analyzer.analyze import analyze
from operon_analyzer.rules import RuleSet, FilterSet


rs = RuleSet().require('transposase') \
              .exclude('cas3') \
              .at_most_n_bp_from_anything('transposase', 500) \
              .same_orientation()

fs = FilterSet().pick_overlapping_features_by_bit_score(0.9) \
                .must_be_within_n_bp_of_anything(1000)

if __name__ == '__main__':
    analyze(sys.stdin, rs, fs)
```

### List of available rules

  * `exclude(feature_name: str)`: Forbid the presence of a particular feature. 
  * `require(feature_name: str)`: Require the presence of a particular feature. 
  * `max_distance(feature1_name: str, feature2_name: str, distance_bp: int)`: The two given features must be no further than `distance_bp` base pairs apart. Requires exactly one of each feature to be present.
  * `at_least_n_bp_from_anything(feature_name: str, distance_bp: int)`: Requires that a feature be at least `distance_bp` base pairs away from any other feature.  This is mostly useful for eliminating overlapping features.
  * `at_most_n_bp_from_anything(feature_name: str, distance_bp: int)`: A given feature must be within `distance_bp` base pairs of another feature. Requires exactly one matching feature to be present. Returns `False` if the given feature is the only feature.
  * `same_orientation(exceptions: Optional[List[str]] = None)`: All features in the operon must have the same orientation.
  * `contains_any_set_of_features(sets: List[List[str]])`: Returns `True` if the operon contains features with all of the names in at least one of the lists. Useful for determining if an operon contains all of the essential genes for a particular system, for example.
  * `contains_exactly_one_of(feature1_name: str, feature2_name: str)`: An exclusive-or of the presence of two features.  That is, one of the features must be present and the other must not.
  * `contains_at_least_n_features(feature_names: List[str], feature_count: int, count_multiple_copies: bool = False)`: The operon must contain at least `feature_count` features in the list. By default, a matching feature that appears multiple times in the operon will only be counted once; to count multiple copies of the same feature, set `count_multiple_copies=True`.
  * `custom(rule: 'Rule')`: Add a rule with a user-defined function. 

### List of available filters

  * `must_be_within_n_bp_of_anything(distance_bp: int)`: If a feature is very far away from anything it's probably not part of an operon.
  * `must_be_within_n_bp_of_feature(feature_name: str, distance_bp: int)`: There may be situations where two features always appear near each other in functional operons.  
  * `pick_overlapping_features_by_bit_score(minimum_overlap_threshold: float)`: If two features overlap by more than `minimum_overlap_threshold`, the one with the lower bit score is ignored.
  * `custom(filt: 'Filter')`: Add a filter with a user-defined function. 

### Analysis Output 

Each line of the CSV will contain an accession ID and the path to the file that contains it, the contig coordinates, and whether it passed or failed the given rules. If it passed, the last column will contain the word `pass` only. Otherwise it will start with `fail` followed by a comma-delimited list of the serialized rules that it failed to adhere to (with the name and parameters that were passed to the method).

## Visualization

Interesting operons can be visualized with a simple gene diagram. It is up to the user to decide how to define this, though this sample script below creates diagrams for all operons that passed all given rules:

```python
import sys
from operon_analyzer.analyze import load_analyzed_operons
from operon_analyzer.visualize import build_operon_dictionary, plot_operons

analysis_csv, pipeline_csv, image_directory = sys.argv[1:]
good_operons = []

with open(pipeline_csv) as f:
    operons = build_operon_dictionary(f)
with open(analysis_csv) as f:
    for contig, filename, start, end, result in load_analyzed_operons(f):
        if result[0] != 'pass':
            continue
        op = operons.get((contig, filename, start, end))
        if op is None:
            continue
        good_operons.append(op)
plot_operons(good_operons, image_directory)
```

## Overview Statistics

Some basic tools are provided to inspect the nature of operons that did not pass all given rules. The intent here is to help researchers determine if their filtering is too aggressive (or not aggressive enough), and to get an overall better feel for the data.

Simple bar plots can be produced as follows:

```python
import sys
import matplotlib.pyplot as plt
from operon_analyzer.analyze import load_analyzed_operons
from operon_analyzer.overview import extract_results, count_results


def plot_bar_chart(filename, title, data, rotate=True):
    fig, ax = plt.subplots()
    x = [str(d[0]).replace(":", "\n") for d in data]
    y = [d[1] for d in data]
    ax.bar(x, y, edgecolor='k')
    if rotate:
        plt.xticks(rotation=90)
    ax.set_ylabel("Count")
    ax.set_title(title)
    plt.savefig("%s.png" % filename, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    lines = load_analyzed_operons(sys.stdin)
    results = extract_results(lines)
    unique_rule_violated, failed_rule_occurrences, rule_failure_counts = count_results(results)
    plot_bar_chart("sole-failure.png", "Number of times that each rule\nwas the only one that failed", sorted(unique_rule_violated.items()))
    plot_bar_chart("total-failures", "Total number of rule failures", sorted(failed_rule_occurrences.items()))
    plot_bar_chart("failures-at-each-contig", "Number of rules failed at each contig", sorted(rule_failure_counts.items()), rotate=False)
```
