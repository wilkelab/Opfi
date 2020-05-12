import csv
from typing import List, Iterator, IO
from operon_analyzer.genes import Operon
from operon_analyzer.rules import RuleSet, Result
from operon_analyzer.parse import parse_pipeline_results, _parse_coordinates


def analyze(input_lines: Iterator[str], ruleset: RuleSet):
    """
    Takes a handle to the CSV from the CRISPR-transposon pipeline
    and user-provided rules, and produces text that describes which
    operons adhered to those rules. If an operon fails any of the rules,
    the exact rules will be enumerated.
    """
    lines = _load_input(input_lines)
    operons = parse_pipeline_results(lines)
    results = _evaluate_operons(operons, ruleset)
    for output_line in _serialize_results(ruleset, results):
        print(output_line)


def load_analyzed_operons(f: IO):
    for line in csv.reader(filter(lambda line: not line.startswith("#"), f)):
        contig, coordinates, result = line
        start, end = _parse_coordinates(coordinates)
        yield contig, start, end, result


def _serialize_results(ruleset: RuleSet, results: List[Result]) -> str:
    """ Generates formatted text for the output of this script. """
    yield "# {rules}".format(rules=str(ruleset))
    for result in results:
        yield str(result)


def _load_input(handle: IO) -> Iterator[List]:
    """ Reads the CSV file produced by the CRISPR-transposon pipeline """
    reader = csv.reader(handle)
    next(reader)  # skip header
    yield from reader


def _evaluate_operons(operons: Iterator[Operon], ruleset: RuleSet) -> Iterator[Result]:
    """ Determines which operons adhere to the filtering rules. """
    for operon in operons:
        yield ruleset.evaluate(operon)
