import csv
from typing import Iterator, IO, Tuple, Optional
from operon_analyzer.genes import Operon
from operon_analyzer.rules import RuleSet, Result, FilterSet
from operon_analyzer.parse import assemble_operons, read_pipeline_output, parse_coordinates


def analyze(input_lines: IO[str], ruleset: RuleSet, filterset: FilterSet = None):
    """
    Takes a handle to the CSV from the CRISPR-transposon pipeline
    and user-provided rules, and produces text that describes which
    operons adhered to those rules. If an operon fails any of the rules,
    the exact rules will be enumerated.
    """
    lines = read_pipeline_output(input_lines)
    operons = assemble_operons(lines)
    results = _evaluate_operons(operons, ruleset, filterset)
    for output_line in _serialize_results(ruleset, results):
        print(output_line)


def load_analyzed_operons(f: IO[str]) -> Iterator[Tuple[str, int, int, str]]:
    """ Loads and parses the data from the output of analyze(). This is
    typically used for analyzing or visualizing candidate operons. """
    for line in csv.reader(filter(lambda line: not line.startswith("#"), f)):
        contig, coordinates, result = line
        start, end = parse_coordinates(coordinates)
        yield contig, start, end, result


def _serialize_results(ruleset: RuleSet, results: Iterator[Result]) -> Iterator[str]:
    """ Generates formatted text for the output of this script. """
    yield "# {rules}".format(rules=str(ruleset))
    for result in results:
        yield str(result)


def _evaluate_operons(operons: Iterator[Operon], ruleset: RuleSet, filterset: Optional[FilterSet] = None) -> Iterator[Result]:
    """ Determines which operons adhere to the filtering rules. """
    for operon in operons:
        if filterset is not None:
            filterset.evaluate(operon)
        yield ruleset.evaluate(operon)
