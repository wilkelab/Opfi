from collections import defaultdict
from typing import Tuple, Dict, DefaultDict, Iterator, IO, List
from operon_analyzer.analyze import load_analyzed_operons


def load_counts(lines: IO[str]) -> Tuple[Dict[str, int],
                                         Dict[str, int],
                                         Dict[int, int]]:
    """ Takes a stream of the output from :func:`operon_analyzer.analyze.analyze` and computes some useful statistics. Namely,
    the number of times each rule was the only one broken for a contig, the number of times each rule was broken regardless of context,
    and the occurrences of the count of broken rules per contig. """
    operons = load_analyzed_operons(lines)
    results = _extract_results(operons)
    unique_rule_violated, failed_rule_occurrences, rule_failure_counts = _count_results(results)
    return unique_rule_violated, failed_rule_occurrences, rule_failure_counts


def _extract_results(lines: Iterator[Tuple[str, int, int, str]]) -> Iterator[str]:
    """ Throws out data we don't care about from the analysis output and returns
    only the actual results """
    for _, _, _, _, result in lines:
        yield result


def _count_results(results: Iterator[List[str]]) -> Tuple[Dict[str, int],
                                                          Dict[str, int],
                                                          Dict[int, int]]:
    """
    Returns a tuple of:
      - the number of times each rule was the only one broken for a contig
      - the number of times each rule was broken regardless of context
      - the occurrences of the count of broken rules per contig
    """

    unique_rule_violated = defaultdict(int)  # type: DefaultDict[str, int]
    failed_rule_occurrences = defaultdict(int)  # type: DefaultDict[str, int]
    rule_failure_counts = defaultdict(int)  # type: DefaultDict[int, int]
    all_rules = set()

    for result in results:
        # result will either be the singular word "pass", or the word "fail" followed
        # by space-separated rules that the operon broke
        if result[0] == 'pass':
            # there were no failed rules
            rule_failure_counts[0] += 1
        else:
            # there was at least one failed rule
            assert result[0] == 'fail'
            for rule in result[1:]:
                failed_rule_occurrences[rule] += 1
                all_rules.add(rule)
            if len(result) == 2:
                # there was exactly one failed rule
                unique_rule_violated[rule] += 1
            # keep track of how many failed rules there were
            rule_failure_counts[len(result[1:])] += 1

    # Ensure rules with zero counts are represented in the data
    for rule in all_rules:
        failed_rule_occurrences[rule] += 0
        unique_rule_violated[rule] += 0
    # convert all results into regular dictionaries since leaving them as defaultdicts
    # could only increase the possibility of a downstream bug
    return dict(unique_rule_violated), dict(failed_rule_occurrences), dict(rule_failure_counts)
