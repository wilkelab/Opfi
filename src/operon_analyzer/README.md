# Operon Analyzer

Tools to filter CRISPR-transposon pipeline results and identify promising candidate operons.

## Analysis

The analysis module provides tools to identify operons that conform to certain rules that interesting candidates adhere to. CSV output is written to stdout, which identifies the outcome of the analysis for each operon.

A sample script that performs this task is given here:

```python
import sys
from operon_analyzer.analyze import analyze
from operon_analyzer.rules import RuleSet


rs = RuleSet().require('transposase') \
              .exclude('cas3') \
              .max_distance_to_anything('transposase', 500) \
              .min_distance_to_anything('transposase', 1)


if __name__ == '__main__':
    analyze(sys.stdin, rs)
```

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
    for contig, start, end, result in load_analyzed_operons(f):
        if result != 'pass':
            continue
        op = operons.get((contig, start, end))
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
