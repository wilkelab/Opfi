from typing import Callable, Optional, List
from operon_analyzer.genes import Feature, Operon


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


class Result(object):
    """
    Records which rules an operon passed and handles serialization of the data.
    Also makes it easy to run follow-up queries.

    """

    def __init__(self, operon: Operon):
        self.operon = operon
        self._passing = []
        self._failing = []

    def add_passing(self, rule: Rule):
        self._passing.append(rule)

    def add_failing(self, rule: Rule):
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

    def min_distance_to_anything(self, feature_name: str, distance_bp: int):
        """
        Requires that a feature be `distance_bp` base pairs from any other feature. 
        This is mostly useful for eliminating overlapping features.
        """
        self._rules.append(Rule('min-distance-to-anything', _min_distance_to_anything, feature_name, distance_bp))
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
        self._rules.append(Rule('same-orientation', _same_orientation, None))
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
    f1 = operon.get(feature1_name)
    f2 = operon.get(feature2_name)
    if len(f1) != 1 or len(f2) != 1:
        return False
    f1 = f1[0]
    f2 = f2[0]
    distance = _feature_distance(f1, f2)
    return 0 <= distance <= distance_bp


def _min_distance_to_anything(operon: Operon, feature_name: str, distance_bp: int) -> bool:
    feature = operon.get(feature_name)
    if len(feature) != 1:
        return False
    feature = feature[0]
    for other_feature in operon:
        if feature is other_feature:
            continue
        distance = _feature_distance(feature, other_feature)
        if distance < distance_bp:
            return False
    return True


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



def _same_orientation(operon: Operon, args=None) -> bool:
    strands = set([feature.strand for feature in operon])
    return len(strands) == 1


def _feature_distance(f1: Feature, f2: Feature) -> int:
    distance1 = f2.start - f1.end
    distance2 = f1.start - f2.end
    # In the case of overlapping features, the distance is defined as 0
    return max(distance1, distance2, 0)
