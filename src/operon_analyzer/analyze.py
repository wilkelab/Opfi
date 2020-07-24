import csv
from collections import defaultdict
from typing import Iterator, IO, Tuple, Optional, List, Dict, Set
from operon_analyzer.genes import Operon
from operon_analyzer.rules import RuleSet, Result, FilterSet
from operon_analyzer.parse import assemble_operons, read_pipeline_output, load_operons
import sys


def analyze(input_lines: IO[str], ruleset: RuleSet, filterset: FilterSet = None, output: IO = None):
    """
    Takes a handle to the CSV from the CRISPR-transposon pipeline
    and user-provided rules, and produces text that describes which
    operons adhered to those rules. If an operon fails any of the rules,
    the exact rules will be enumerated.
    """
    output = sys.stdout if output is None else output
    lines = read_pipeline_output(input_lines)
    operons = assemble_operons(lines)
    results = _evaluate_operons(operons, ruleset, filterset)
    output.write("# {rules}\n".format(rules=str(ruleset)))
    writer = csv.writer(output)
    for result in results:
        line = [result.operon.contig, result.operon.contig_filename, result.operon.start, result.operon.end]
        if result.is_passing:
            line.append("pass")
        else:
            line.append("fail")
            for rule in result._failing:
                line.append(str(rule))
        writer.writerow(line)


def load_analyzed_operons(f: IO[str]) -> Iterator[Tuple[str, int, int, str]]:
    """ Loads and parses the data from the output of analyze(). This is
    typically used for analyzing or visualizing candidate operons. """
    for line in csv.reader(filter(lambda line: not line.startswith("#"), f)):
        contig, contig_filename, start, end = line[:4]
        start = int(start)
        end = int(end)
        result = line[4:]
        yield contig, contig_filename, start, end, result


def _evaluate_operons(operons: Iterator[Operon], ruleset: RuleSet, filterset: Optional[FilterSet] = None) -> Iterator[Result]:
    """ Determines which operons adhere to the filtering rules. """
    for operon in operons:
        if filterset is not None:
            filterset.evaluate(operon)
        yield ruleset.evaluate(operon)


def cluster_operons_by_feature_order(operons: Iterator[Operon]):
    """ Organizes all operons into a dictionary based on the order/identity of the features.
    Cases where the overall order is inverted are considered to be the same. The keys of the dictionary
    are the dash-delimited feature names, with one of the two orientations (if both exist) arbitrarily chosen.
    If there are ignored features, they will not appear in the key. """
    bins = defaultdict(list)
    for operon in operons:
        feature_names = _get_sorted_feature_names(operon)
        reverse_feature_names = tuple(reversed(feature_names))
        if reverse_feature_names not in bins:
            bins[feature_names].append(operon)
        else:
            bins[reverse_feature_names].append(operon)
    return bins


def _get_diffed_cluster_keys(clustered_operons: Dict[str, Operon], diff_against: Dict[str, Operon]) -> Set[str]:
    """ Given two sets of clustered operons, we want to know what the unique keys are that exist in `clustered_operons` so we can later only plot those. """
    return set(clustered_operons.keys()) - set(diff_against.keys())


def _get_sorted_feature_names(operon: Operon) -> List[str]:
    """ Produces a list of feature names, ordered by the lowest genomic coordinate
    (that is, regardless of the orientation of the feature, whichever end has the lowest
    numbered nucleotide position is what is used here). """
    return tuple((feature.name for feature in sorted(operon, key=lambda feat: feat.start)))
