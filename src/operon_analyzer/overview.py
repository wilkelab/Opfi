#    plot number of failed rules
#    for those with only one rule broken, counts of each rule violated
#    counts of all rules violated

from collections import defaultdict
import sys
from typing import Tuple, DefaultDict, Iterator
import matplotlib.pyplot as plt
from operon_analyzer.analyze import load_analyzed_operons


def extract_results(lines: Iterator[Tuple[str, int, int, str]]) -> Iterator[str]:
    for _, _, _, result in lines:
        yield result


def count_results(results: Iterator[str]) -> Tuple[DefaultDict[str, int],
                                                   DefaultDict[str, int],
                                                   DefaultDict[str, int]]:
    """
    Returns a tuple of:
      - the number of times each rule was the only one broken for a contig
      - the number of times each rule was broken regardless of context
      - the occurrences of the count of broken rules per contig
    """

    unique_rule_violated = defaultdict(int)
    failed_rule_occurrences = defaultdict(int)
    rule_failure_counts = defaultdict(int)
    all_rules = set()

    for result in results:
        data = result.split()
        if data[0] == 'pass':
            rule_failure_counts[0] += 1
        else:
            assert data[0] == 'fail'
            for rule in data[1:]:
                failed_rule_occurrences[rule] += 1
                all_rules.add(rule)
            if len(data) == 2:
                unique_rule_violated[rule] += 1
            rule_failure_counts[len(data[1:])] += 1
    for rule in all_rules:
        # Ensure rules with zero counts are represented in the data
        failed_rule_occurrences[rule] += 0
        unique_rule_violated[rule] += 0
    return unique_rule_violated, failed_rule_occurrences, rule_failure_counts


if __name__ == '__main__':
    lines = load_analyzed_operons(sys.stdin)
    results = extract_results(lines)
    unique_rule_violated, failed_rule_occurrences, rule_failure_counts = count_results(results)
    print(unique_rule_violated)
    print(failed_rule_occurrences)
    print(rule_failure_counts)
