import csv
from collections import defaultdict
from typing import List, Iterator, IO, Tuple, Callable, Union, Any, Optional


def analyze(input_lines: Iterator[str], ruleset: 'RuleSet'):
    """
    Takes a handle to the CSV from the CRISPR-transposon pipeline
    and user-provided rules, and produces text that describes which
    operons adhered to those rules. If an operon fails any of the rules,
    the exact rules will be enumerated.
    """
    lines = _load_input(input_lines)
    operons = _parse_pipeline_results(lines)
    results = _evaluate_operons(operons, ruleset)
    for output_line in _serialize_results(ruleset, results):
        print(output_line)


class Feature(object):
    """ Represents a gene or CRISPR repeat array. """

    def __init__(self,
                 name: str,
                 coordinates: Tuple[int, int],
                 orfid: str,
                 strand: int,
                 accession: str,
                 e_val: float,
                 description: str,
                 sequence: str):
        self.name = name
        self.coordinates = coordinates
        self.orfid = orfid
        self.strand = strand
        self.accession = accession
        self.e_val = e_val
        self.description = description
        self.sequence = sequence

    @property
    def start(self) -> int:
        """
        The leftmost position of the Feature, relative to the
        orientation of the Operon.
        """
        return min(self.coordinates)

    @property
    def end(self) -> int:
        """
        The rightmost position of the Feature, relative to the
        orientation of the Operon.
        """
        return max(self.coordinates)


class Operon(object):
    """
    Provides access to Features that were found in the same genomic region,
    which presumably comprise an actual operon. Whether this is true in reality
    must be determined by the user, if that is meaningful to them.
    """

    def __init__(self,
                 contig: str,
                 start: int,
                 end: int,
                 features: List[Feature]):
        assert contig, "Missing contig name"
        assert 0 <= start and 0 <= end, "Invalid contig position"
        assert features, "Contig did not contain any features"
        self.contig = contig
        self.start = start
        self.end = end
        self._features = features

    def __iter__(self):
        yield from self._features

    @property
    def feature_names(self):
        """ Iterates over the name of each feature in the operon """
        yield from (feature.name for feature in self._features)

    def get(self, feature_name: str) -> List[Feature]:
        """ Returns a list of every Feature with a given name. """
        features = []
        for feature in self._features:
            if feature.name == feature_name:
                features.append(feature)
        return features


class Result(object):
    """
    Records which rules an operon passed and handles serialization of the data.
    Also makes it easy to run follow-up queries.

    """

    def __init__(self, operon: Operon):
        self.operon = operon
        self._passing = []
        self._failing = []

    def add_passing(self, rule: 'Rule'):
        self._passing.append(rule)

    def add_failing(self, rule: 'Rule'):
        self._failing.append(rule)

    @property
    def is_passing(self) -> bool:
        """ Declares whether the given Operon adhered to all given Rules. """
        assert (self._passing or self._failing), "Result has no Rules!"
        return bool(self._passing) and not bool(self._failing)

    def __repr__(self) -> str:
        outcome = "pass" if self.is_passing else "fail %s" % ",".join(str(rule) for rule in self._failing)
        text = "{contig},{start}..{end},{outcome}"
        return text.format(
                contig=self.operon.contig,
                start=self.operon.start,
                end=self.operon.end,
                outcome=outcome)


class Rule(object):
    """ Defines a requirement that elements of an operon must adhere to. """

    def __init__(self, name: str, function: Callable, *args):
        self._name = name
        self._function = function
        self._args = args

    def evaluate(self, operon: Operon) -> bool:
        return self._function(operon, *self._args)

    def __repr__(self) -> str:
        return "{name}:{args}".format(
                name=str(self._name),
                args="-".join(map(str, self._args)))


class RuleSet(object):
    """ Creates, stores and evaluates `Rule`s that an operon must adhere to."""

    def __init__(self):
        self._rules = []

    def exclude(self, feature_name: str):
        """ Forbid the presence of a particular feature. """
        def func(operon: Operon, feature_name: str) -> bool:
            return feature_name not in operon.feature_names
        self._rules.append(Rule('exclude', func, feature_name))
        return self

    def require(self, feature_name: str):
        """ Require the presence of a particular feature. """
        def func(operon: Operon, feature_name: str) -> bool:
            return feature_name in operon.feature_names
        self._rules.append(Rule('require', func, feature_name))
        return self

    def max_distance(self, feature1_name: str, feature2_name: str, distance_bp: int):
        """
        The two given features must be no further than distance_bp base pairs
        apart. Requires exactly one of each feature to be present.
        """

        def func(operon: Operon, feature1_name: str, feature2_name: str, distance_bp: int) -> bool:
            f1 = operon.get(feature1_name)
            f2 = operon.get(feature2_name)
            if len(f1) != 1 or len(f2) != 1:
                return False
            f1 = f1[0]
            f2 = f2[0]
            distance = _feature_distance(f1, f2)
            return 0 <= distance <= distance_bp
        self._rules.append(Rule('max-distance', func, feature1_name, feature2_name, distance_bp))
        return self

    def max_distance_to_anything(self, feature_name: str, distance_bp: int):
        """
        A given feature must be within distance_bp base pairs of any other feature.
        Requires exactly one matching feature to be present.
        Returns False if the given feature is the only feature.
        """
        self._rules.append(Rule('max-distance-to-anything', _max_distance_to_anything, feature_name, distance_bp))
        return self

    def same_orientation(self, exceptions: Optional[List[str]] = None):
        """
        All features in the operon must have the same orientation.
        """
        def func(operon: Operon, args=None) -> bool:
            strands = set([feature.strand for feature in operon])
            return len(strands) == 1
        self._rules.append(Rule('same-orientation', func, None))
        return self

    def custom(self, rule: 'Rule'):
        """ Add a rule with a user-defined function. """
        self._rules.append(rule)
        return self

    def evaluate(self, operon: Operon) -> Result:
        """ See if an operon adheres to all rules. """
        assert self._rules, "RuleSet does not have any Rules"
        result = Result(operon)
        for rule in self._rules:
            if rule.evaluate(operon):
                result.add_passing(rule)
            else:
                result.add_failing(rule)
        return result

    def __repr__(self) -> str:
        return ",".join((str(rule) for rule in self._rules))


def _max_distance_to_anything(operon: Operon, feature_name: str, distance_bp: int) -> bool:
    feature = operon.get(feature_name)
    if len(feature) != 1:
        return False
    feature = feature[0]
    for other_feature in operon:
        if feature is other_feature:
            continue
        distance = _feature_distance(feature, other_feature) 
        if distance <= distance_bp:
            return True
    return False


def _feature_distance(f1: Feature, f2: Feature) -> int:
    distance1 = f2.start - f1.end
    distance2 = f1.start - f2.end
    assert (distance1 <= 0 or distance2 <= 0) and (distance1 >= 0 or distance2 >= 0)
    return max(distance1, distance2)


def _parse_coordinates(coordinates: str) -> Tuple[int, int]:
    """ Parses base pair ranges denoting position in a contig. """
    start, end = coordinates.split("..")
    return int(start), int(end)


def _parse_feature(line: List[str]) -> (str, Tuple[int, int], Feature):
    """ Creates a Feature from a line of output from a CSVReader """
    contig = line[0]
    coordinates = _parse_coordinates(line[1])
    feature = line[2]
    feature_coordinates = _parse_coordinates(line[3])
    query_orfid = line[4]
    strand = int(line[5]) if line[5] else None
    hit_accession = line[6]
    hit_eval = float(line[7]) if line[7] else None
    description, sequence = line[8:]
    return contig, coordinates, Feature(
        feature,
        feature_coordinates,
        query_orfid,
        strand,
        hit_accession,
        hit_eval,
        description,
        sequence)


def _parse_pipeline_results(lines: Iterator[Tuple]) -> Iterator['Operon']:
    """
    Takes the output from the CRISPR-transposon pipeline, loads all features,
    and assembles them into putative Operons.
    """
    operon_features = defaultdict(list)
    for line in lines:
        contig, coordinates, feature = _parse_feature(line)
        operon_features[(contig, coordinates)].append(feature)
    for (contig, (start, end)), features in operon_features.items():
        yield Operon(contig, start, end, features)


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
