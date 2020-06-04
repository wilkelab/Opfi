from typing import Callable, Optional, List
from operon_analyzer.genes import Feature, Operon


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

    def __repr__(self) -> str:
        outcome = "pass" if self.is_passing else "fail %s" % " ".join(str(rule) for rule in self._failing)
        text = "{contig},{start}..{end},{outcome}"
        return text.format(
                contig=self.operon.contig,
                start=self.operon.start,
                end=self.operon.end,
                outcome=outcome)


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

    def must_be_within_n_bp_of_feature(self, feature_name: str, distance_bp: int):
        """ There may be situations where two features always appear near each other in functional operons. """
        self._filters.append(Filter('must-be-within-n-bp-of-feature', _must_be_within_n_bp_of_feature, feature_name, distance_bp))
        return self

    def pick_overlapping_features_by_bit_score(self, minimum_overlap_threshold: float):
        """ If two features overlap by more than `minimum_overlap_threshold`, the one with the lower bit score is ignored. """
        self._filters.append(Filter('overlaps-other-feature', _pick_overlapping_features_by_bit_score, minimum_overlap_threshold))
        return self

    def evaluate(self, operon: Operon):
        """ Run the filters on the operon and set Features that fail to meet the requirements to be ignored. """
        for filt in self._filters:
            filt.run(operon)


def _pick_overlapping_features_by_bit_score(operon: Operon, minimum_overlap_threshold: float):
    # build lists of overlapping feature groups
    # for a case where A overlaps with B, and B overlaps with C, but A and C don't overlap, what should happen?
    # for each feature, make a list of other features that overlap by at least N%
    # iterate over that list ignore whichever feature is worse, unless the other feature has already been ignored.
    # This last bit is necessary since we don't want a situation where all our features just disappear
    pass


def _must_be_within_n_bp_of_anything(operon: Operon, ignored_reason_message: str, distance_bp: int):
    for feature in operon.all_features:
        distances = _calculate_all_distances(operon, feature.name)
        if min(distances) > distance_bp:
            feature.ignore(ignored_reason_message)


def _must_be_within_n_bp_of_feature(operon: Operon, ignored_reason_message: str, feature_name: str, distance_bp: int):
    for feature in operon.all_features:
        if feature.name == feature_name:
            continue
        if not _max_distance(operon, feature_name, feature.name, distance_bp):
            feature.ignore(ignored_reason_message)


class RuleSet(object):
    """ Creates, stores and evaluates `Rule`s that an operon must adhere to."""

    def __init__(self):
        self._rules = []

    def exclude(self, feature_name: str):
        """ Forbid the presence of a particular feature. """
        self._rules.append(Rule('exclude', _exclude, feature_name))
        return self

    def require(self, feature_name: str):
        """ Require the presence of a particular feature. """
        self._rules.append(Rule('require', _require, feature_name))
        return self

    def max_distance(self, feature1_name: str, feature2_name: str, distance_bp: int):
        """
        The two given features must be no further than distance_bp base pairs
        apart. Requires exactly one of each feature to be present.
        """
        self._rules.append(Rule('max-distance', _max_distance, feature1_name, feature2_name, distance_bp))
        return self

    def at_least_n_bp_from_anything(self, feature_name: str, distance_bp: int):
        """
        Requires that a feature be at least `distance_bp` base pairs away from any other feature.
        This is mostly useful for eliminating overlapping features.
        """
        self._rules.append(Rule('at-least-n-bp-from-anything', _at_least_n_bp_from_anything, feature_name, distance_bp))
        return self

    def at_most_n_bp_from_anything(self, feature_name: str, distance_bp: int):
        """
        A given feature must be within distance_bp base pairs of another feature.
        Requires exactly one matching feature to be present.
        Returns False if the given feature is the only feature.
        """
        self._rules.append(Rule('at-most-n-bp-from-anything', _at_most_n_bp_from_anything, feature_name, distance_bp))
        return self

    def same_orientation(self, exceptions: Optional[List[str]] = None):
        """
        All features in the operon must have the same orientation.
        """
        self._rules.append(Rule('same-orientation', _same_orientation, None))
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

    def contains_exactly_one_of(self, feature1_name: str, feature2_name: str):
        """
        An exclusive-or of the presence of two features.
        That is, one of the features must be present and the other must not.
        """
        self._rules.append(Rule('contains-exactly-one-of',
                                _contains_exactly_one_of,
                                feature1_name,
                                feature2_name))
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


def _exclude(operon: Operon, feature_name: str) -> bool:
    """ Returns false if a feature's name in the operon matches the given string. Case insensitive. """
    return feature_name not in operon.feature_names


def _require(operon: Operon, feature_name: str) -> bool:
    """ Returns true if a feature's name in the operon matches the given string. Case insensitive. """
    return feature_name in map(str.lower, operon.feature_names)


def _max_distance(operon: Operon, feature1_name: str, feature2_name: str, distance_bp: int) -> bool:
    """ Returns whether two given Features are within distance_bp base pairs from each other. """
    f1 = operon.get_unique(feature1_name)
    f2 = operon.get_unique(feature2_name)
    if f1 is None or f2 is None:
        return False
    distance = _feature_distance(f1, f2)
    return 0 <= distance <= distance_bp


def _calculate_all_distances(operon: Operon, feature_name: str) -> Optional[int]:
    """ Calculates the distance of a unique feature to all other features. """
    distances = []
    feature = operon.get_unique(feature_name)
    if feature is None:
        return None
    for other_feature in operon:
        if feature is other_feature:
            continue
        distance = _feature_distance(feature, other_feature)
        distances.append(distance)
    return distances


def _at_least_n_bp_from_anything(operon: Operon, feature_name: str, distance_bp: int) -> bool:
    """ Whether a given feature is more than distance_bp base pairs from another Feature. """
    distances = _calculate_all_distances(operon, feature_name)
    if distances is None:
        return False
    if not distances:
        return True
    return min(distances) >= distance_bp


def _at_most_n_bp_from_anything(operon: Operon, feature_name: str, distance_bp: int) -> bool:
    """ Whether a given feature is less than distance_bp base pairs from any other feature. """
    distances = _calculate_all_distances(operon, feature_name)
    if distances is None:
        return False
    return min(distances) <= distance_bp


def _same_orientation(operon: Operon, args=None) -> bool:
    """ Whether every gene is transcribed in the same direction. """
    strands = set([feature.strand for feature in operon])
    return len(strands) == 1


def _contains_features(operon: Operon, feature_names: List[str]) -> bool:
    """ Whether an Operon contains features with the given names. """
    return len(set(operon.feature_names) & set(feature_names)) == len(feature_names)


def _contains_any_set_of_features(operon: Operon, sets: List[List[str]]) -> bool:
    """ Whether any the Operon has features with any of the given sets of names. """
    return any([_contains_features(operon, feature_names) for feature_names in sets])


def _contains_exactly_one_of(operon: Operon, f1: str, f2: str) -> bool:
    """ Whether the operon has one feature or another, but not both. """
    return (f1 in operon.feature_names) ^ (f2 in operon.feature_names)


def _feature_distance(f1: Feature, f2: Feature) -> int:
    """ Returns the distance between two Features in base pairs. """
    distance1 = f2.start - f1.end
    distance2 = f1.start - f2.end
    # In the case of overlapping features, the distance is defined as 0
    return max(distance1, distance2, 0)
