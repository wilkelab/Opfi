import re
from typing import Callable, Optional, List
from operon_analyzer.genes import Feature, Operon
import more_itertools


class SerializableFunction(object):
    """ A base class for functions that we need to be able to serialize. Do not instantiate this directly. """

    def __init__(self, name: str, function: Callable, *args, custom_repr: Optional[str] = None):
        self._name = name
        self._function = function
        self._args = args
        self._custom_repr = custom_repr

    def __repr__(self) -> str:
        if self._custom_repr is not None:
            return self._custom_repr
        return "{name}:{args}".format(
                name=str(self._name),
                args="-".join(map(str, self._args)))


class Rule(SerializableFunction):
    """ Defines a requirement that elements of an operon must adhere to. """

    def evaluate(self, operon: Operon) -> bool:
        """ Determine if an operon adheres to this rule. """
        return self._function(operon, *self._args)


class Result(object):
    """
    Records which rules an operon passed and handles serialization of the data.
    Also makes it easy to run follow-up queries.

    """

    def __init__(self, operon: Operon):
        self.operon = operon
        self._passing = []  # type: List[Rule]
        self._failing = []  # type: List[Rule]

    def add_passing(self, rule: Rule):
        """ Mark this rule as being one that the operon passed. """
        self._passing.append(rule)

    def add_failing(self, rule: Rule):
        """ Mark this rule as being one that the operon failed. """
        self._failing.append(rule)

    @property
    def is_passing(self) -> bool:
        """ Declares whether the given Operon adhered to all given Rules. """
        assert (self._passing or self._failing), "Result has no Rules!"
        return bool(self._passing) and not bool(self._failing)


class Filter(SerializableFunction):
    """
    A function that will be run on an Operon that marks Features as being ignorable
    for the purposes of evaluating RuleSets.
    """

    def run(self, operon: Operon):
        """ Mark Features in the Operon as ignored if they don't pass the filter. This
        will prevent them from being taken into account during Rule evaluation and (by
        default) during visualization. """
        return self._function(operon, str(self), *self._args)


class FilterSet(object):
    """
    Stores functions that take an Operon and mark individual Features as ignored
    in case we think they are not actually worth taking into account when evaluating rules.
    Features can be ignored for multiple reasons.
    """

    def __init__(self):
        self._filters = []

    def must_be_within_n_bp_of_anything(self, distance_bp: int):
        """ If a feature is very far away from anything it's probably not part of an operon. """
        self._filters.append(Filter('must-be-within-n-bp-of-anything', _must_be_within_n_bp_of_anything, distance_bp))
        return self

    def must_be_within_n_bp_of_feature(self, feature_name: str, distance_bp: int, regex: bool = False):
        """ There may be situations where two features always appear near each other in functional operons. """
        self._filters.append(Filter('must-be-within-n-bp-of-feature', _must_be_within_n_bp_of_feature, feature_name, distance_bp, regex))
        return self

    def pick_overlapping_features_by_bit_score(self, minimum_overlap_threshold: float):
        """ If two features overlap by more than ``minimum_overlap_threshold``, the one with the lower bit score is ignored. """
        self._filters.append(Filter('overlaps-%s',
                                    _pick_overlapping_features_by_bit_score,
                                    minimum_overlap_threshold))
        return self

    def custom(self, filt: 'Filter'):
        """ Add a rule with a user-defined function. """
        self._filters.append(filt)
        return self

    def evaluate(self, operon: Operon):
        """ 
        Run the filters on the operon and set Features that fail to meet the requirements to be ignored. 

        Args:
            operon: The :class:`operon_analyzer.genes.Operon` object whose features will be evaluated.
        """
        for filt in self._filters:
            filt.run(operon)


def _pick_overlapping_features_by_bit_score(operon: Operon, ignored_reason_message: str, minimum_overlap_threshold: float):
    """ If BLAST identified two genes at the same location, we want to determine which one we should
    actually consider the ORF to actually be. In these cases, we pick whichever one has the highest
    bit score. Since there are many small overlaps, the user must specify some arbitrary limit to
    define what an overlap is. """

    for feature in operon.all_genes:
        for other_feature in operon.all_genes:
            if feature is other_feature:
                # don't compare feature to itself
                continue
            overlap = _calculate_overlap(feature, other_feature)
            if overlap is None:
                # these features do not overlap
                continue
            if overlap >= minimum_overlap_threshold and feature.bit_score < other_feature.bit_score:
                feature.ignore(ignored_reason_message % other_feature.name)


def _calculate_overlap(feature: Feature, other_feature: Feature) -> Optional[float]:
    """ Calculates the fraction of a feature that overlaps with another feature. """
    if feature.start > other_feature.end or other_feature.start > feature.end:
        # these features don't overlap at all
        return None
    # Find the lower of the two ends
    end = min(feature.end, other_feature.end)
    # Find the higher of the two starts
    start = max(feature.start, other_feature.start)
    # Determine how much overlap there is
    feature_length = feature.end - feature.start + 1
    return (end - start + 1) / feature_length


def _must_be_within_n_bp_of_anything(operon: Operon, ignored_reason_message: str, distance_bp: int):
    for feature in operon.all_features:
        distances = _calculate_all_distances(operon, feature)
        if not distances:
            continue
        if min(distances) > distance_bp:
            feature.ignore(ignored_reason_message)


def _must_be_within_n_bp_of_feature(operon: Operon, ignored_reason_message: str, feature_name: str, distance_bp: int, regex: bool):
    reg = re.compile(feature_name, re.IGNORECASE)
    for feature in operon.all_features:
        if reg.match(feature.name):
            continue
        if not _max_distance(operon, feature_name, feature.name, distance_bp, False, regex):
            feature.ignore(ignored_reason_message)


class RuleSet(object):
    """ Creates, stores and evaluates :class:`operon_analyzer.rules.Rule` s that an operon must adhere to."""

    def __init__(self):
        self._rules = []

    def exclude(self, feature_name: str, regex: bool = False):
        """ Forbid the presence of a particular feature. """
        self._rules.append(Rule('exclude', _exclude, feature_name, regex))
        return self

    def require(self, feature_name: str, regex: bool = False):
        """ Require the presence of a particular feature. """
        self._rules.append(Rule('require', _require, feature_name, regex))
        return self

    def max_distance(self, feature1_name: str, feature2_name: str, distance_bp: int, closest_pair_only: bool = False, regex: bool = False):
        """
        The two given features must be no further than ``distance_bp`` base pairs
        apart. If there is more than one match, all possible pairs must meet the criteria,
        unless ``closest_pair_only`` is True in which case only the closets pair is considered.
        """
        self._rules.append(Rule('max-distance', _max_distance, feature1_name, feature2_name, distance_bp, closest_pair_only, regex))
        return self

    def at_least_n_bp_from_anything(self, feature_name: str, distance_bp: int, regex=False):
        """
        Requires that a feature be at least ``distance_bp`` base pairs away from any other feature.
        This is mostly useful for eliminating overlapping features.
        """
        self._rules.append(Rule('at-least-n-bp-from-anything', _at_least_n_bp_from_anything, feature_name, distance_bp, regex))
        return self

    def at_most_n_bp_from_anything(self, feature_name: str, distance_bp: int, regex: bool = False):
        """
        A given feature must be within ``distance_bp`` base pairs of another feature.
        Requires exactly one matching feature to be present.
        Returns False if the given feature is the only feature.
        """
        self._rules.append(Rule('at-most-n-bp-from-anything', _at_most_n_bp_from_anything, feature_name, distance_bp, regex))
        return self

    def same_orientation(self, exceptions: Optional[List[str]] = None):
        """
        All features in the operon must have the same orientation.
        """
        self._rules.append(Rule('same-orientation', _same_orientation, exceptions))
        return self

    def contains_any_set_of_features(self, sets: List[List[str]]):
        """
        Returns True if the operon contains features with all of the names
        in at least one of the lists. Useful for determining if an operon contains
        all of the essential genes for a particular system, for example.
        """
        serialized_sets = "|".join(["-".join(names) for names in sets])
        custom_repr = f'contains-any-set-of-features:{serialized_sets}'
        self._rules.append(Rule('contains-any-set-of-features',
                                _contains_any_set_of_features,
                                sets, custom_repr=custom_repr))
        return self

    def contains_exactly_one_of(self, feature1_name: str, feature2_name: str, regex: bool = False):
        """
        An exclusive-or of the presence of two features.
        That is, one of the features must be present and the other must not.
        """
        self._rules.append(Rule('contains-exactly-one-of',
                                _contains_exactly_one_of,
                                feature1_name,
                                feature2_name))
        return self

    def contains_at_least_n_features(self, feature_names: List[str], feature_count: int, count_multiple_copies: bool = False):
        """
        The operon must contain at least ``feature_count`` features in the list. By default, a
        matching feature that appears multiple times in the operon will only be counted once;
        to count multiple copies of the same feature, set ``count_multiple_copies`` to True.
        """
        serialized_list = "|".join(feature_names)
        custom_repr = f'contains-at-least-n-features:{serialized_list}-{feature_count}-{count_multiple_copies}'
        self._rules.append(Rule('contains_at_least_n_features',
                                _contains_at_least_n_features,
                                feature_names,
                                feature_count,
                                count_multiple_copies,
                                custom_repr=custom_repr))
        return self

    def contains_group(self, feature_names: List[str], max_gap_distance_bp: int, require_same_orientation: bool):
        """ The operon must contain a contiguous set of features (in any order) separated by no more than ``max_gap_distance_bp``.
        Optionally, the user may require that the features must all have the same orientation. """
        self._rules.append(Rule('contains_group',
                                _contains_group,
                                feature_names,
                                max_gap_distance_bp,
                                require_same_orientation))
        return self

    def maximum_size(self, feature_name: str, max_bp: int, all_matching_features_must_pass: bool = False, regex: bool = False):
        """ The operon must contain at least one feature with ``feature_name`` with a size (in base pairs) of ``max_bp`` or smaller.
        If ``all_matching_features_must_pass`` is True, every matching Feature must be at least ``max_bp`` long. """
        self._rules.append(Rule('maximum_size',
                                _maximum_size,
                                feature_name,
                                max_bp,
                                all_matching_features_must_pass,
                                regex))
        return self

    def minimum_size(self, feature_name: str, min_bp: int, all_matching_features_must_pass: bool = False, regex: bool = False):
        """ The operon must contain at least one feature with ``feature_name`` with a size (in base pairs) of ``min_bp`` or larger.
        If ``all_matching_features_must_pass`` is True, every matching Feature must be at least ``min_bp`` long. """
        self._rules.append(Rule('minimum_size',
                                _minimum_size,
                                feature_name,
                                min_bp,
                                all_matching_features_must_pass,
                                regex))
        return self

    def custom(self, rule: 'Rule'):
        """ Add a rule with a user-defined function. """
        self._rules.append(rule)
        return self

    def evaluate(self, operon: Operon) -> Result:
        """ 
        See if an operon adheres to all rules.

        Args:
            operon: The :class:`operon_analyzer.genes.Operon` object to evaluate.
        """
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


def _exclude(operon: Operon, feature_name: str, regex: str) -> bool:
    """ Returns false if a feature's name in the operon matches the given string. """
    return not _require(operon, feature_name, regex)


def _require(operon: Operon, feature_name: str, regex: bool) -> bool:
    """ Returns true if a feature's name in the operon matches the given string. """
    if regex:
        rx = re.compile(feature_name, re.IGNORECASE)
        return any([rx.match(name) for name in operon.feature_names])
    return feature_name.lower() in map(str.lower, operon.feature_names)


def _max_distance(operon: Operon, feature1_name: str, feature2_name: str, distance_bp: int, closest_pair_only: bool, regex: bool) -> bool:
    """ Returns whether two given Features are within distance_bp base pairs from each other.
    This must hold for all copies of the features match unless closest_pair_only is True,
    in which case we only consider the pair that is the shortest distance apart. """
    distances = []
    for f1 in operon.get(feature1_name, regex):
        for f2 in operon.get(feature2_name, regex):
            if f1 is f2:
                continue
            distance = _feature_distance(f1, f2)
            distances.append(distance)
    if not distances:
        return False
    if closest_pair_only:
        return min(distances) <= distance_bp
    return all([distance <= distance_bp for distance in distances])


def _calculate_all_distances(operon: Operon, feature: Feature) -> Optional[int]:
    """ Calculates the distance of a unique feature to all other features. """
    distances = []
    for other_feature in operon:
        if feature is other_feature:
            continue
        distance = _feature_distance(feature, other_feature)
        distances.append(distance)
    return distances


def _at_least_n_bp_from_anything(operon: Operon, feature_name: str, distance_bp: int, regex: bool = False) -> bool:
    """ Whether a given feature is more than distance_bp base pairs from another Feature. """
    at_least_one_good = False
    for feature in operon.get(feature_name, regex):
        distances = _calculate_all_distances(operon, feature)
        if not distances:
            return True
        if min(distances) < distance_bp:
            return False
        at_least_one_good = True
    return at_least_one_good


def _at_most_n_bp_from_anything(operon: Operon, feature_name: str, distance_bp: int, regex: bool = False) -> bool:
    """ Whether a given feature is less than distance_bp base pairs from any other feature. """
    at_least_one_good = False
    for feature in operon.get(feature_name, regex):
        distances = _calculate_all_distances(operon, feature)
        if not distances:
            return False
        if min(distances) > distance_bp:
            return False
        at_least_one_good = True
    return at_least_one_good


def _same_orientation(operon: Operon, exceptions: Optional[List[str]]) -> bool:
    """ Whether every gene is transcribed in the same direction. """
    exceptions = [] if exceptions is None else exceptions
    strands = set([feature.strand for feature in operon if feature.name.lower() not in map(str.lower, exceptions)])
    return len(strands) == 1


def _contains_features(operon: Operon, feature_names: List[str]) -> bool:
    """ Whether an Operon contains features with the given names. """
    return len(set(map(str.lower, operon.feature_names)) & set(map(str.lower, feature_names))) == len(feature_names)


def _contains_any_set_of_features(operon: Operon, sets: List[List[str]]) -> bool:
    """ Whether any the Operon has features with any of the given sets of names. """
    return any([_contains_features(operon, feature_names) for feature_names in sets])


def _contains_exactly_one_of(operon: Operon, f1: str, f2: str) -> bool:
    """ Whether the operon has one feature or another, but not both. """
    reg1 = re.compile(f1, re.IGNORECASE)
    reg2 = re.compile(f2, re.IGNORECASE)
    return any(map(reg1.match, operon.feature_names)) ^ any(map(reg2.match, operon.feature_names))


def _contains_at_least_n_features(operon: Operon, feature_names: List[str], feature_count: int, count_multiple_copies: bool) -> bool:
    """ Whether the operon has at least feature_count given features. """
    matches = [feature_name for feature_name in operon.feature_names if feature_name.lower() in map(str.lower, feature_names)]
    if len(matches) >= feature_count and count_multiple_copies:
        return True
    elif len(set(matches)) >= feature_count:
        return True
    else:
        return False


def _feature_distance(f1: Feature, f2: Feature) -> int:
    """ Returns the distance between two Features in base pairs. """
    distance1 = f2.start - f1.end
    distance2 = f1.start - f2.end
    # In the case of overlapping features, the distance is defined as 0
    return max(distance1, distance2, 0)


def _contains_group(operon: Operon, feature_names: List[str], max_gap_distance_bp: int, require_same_orientation: bool) -> bool:
    """ Determines whether Features with exactly the names in feature_names occur in any order without interruption by other Features, with gaps between each Feature no larger than max_gap_distance_bp. Optionally, every Feature can be required to be in the same orientation. """
    assert len(feature_names) > 1
    if len(operon) < len(feature_names):
        return False

    sorted_feature_names = tuple(sorted(map(str.lower, feature_names)))
    for operon_chunk in more_itertools.windowed(operon, len(feature_names)):
        sorted_chunk_names = tuple(sorted([feature.name.lower() for feature in operon_chunk]))
        if sorted_chunk_names == sorted_feature_names:
            max_gap_distance_in_group = 0
            for feature1, feature2 in zip(operon_chunk, operon_chunk[1:]):
                max_gap_distance_in_group = max(_feature_distance(feature1, feature2), max_gap_distance_in_group)
            if max_gap_distance_in_group > max_gap_distance_bp:
                continue
            if require_same_orientation:
                strands = set((feature.strand for feature in operon_chunk if feature.strand is not None))
                if len(strands) != 1:
                    continue
            return True
    return False


def _minimum_size(operon: Operon, feature_name: str, min_bp: int, all_matching_features_must_pass: bool, regex: bool) -> bool:
    results = [len(feature) >= min_bp for feature in operon.get(feature_name, regex)]
    if all_matching_features_must_pass:
        return all(results)
    return any(results)


def _maximum_size(operon: Operon, feature_name: str, max_bp: int, all_matching_features_must_pass: bool, regex: bool) -> bool:
    results = [len(feature) <= max_bp for feature in operon.get(feature_name, regex)]
    if all_matching_features_must_pass:
        return all(results)
    return any(results)
