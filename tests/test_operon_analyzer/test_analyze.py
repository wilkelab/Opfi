from operon_analyzer.genes import Feature, Operon
from operon_analyzer.rules import RuleSet, _feature_distance, _max_distance, _contains_features
from operon_analyzer.analyze import _serialize_results
from operon_analyzer.visualize import calculate_adjusted_operon_bounds
from operon_analyzer.overview import _count_results
import pytest
from hypothesis.strategies import composite, text, integers, sampled_from, floats, lists
from hypothesis import given, settings
import string

name_characters = string.ascii_lowercase + string.ascii_uppercase + string.digits

sequence_characters = 'ACDEFGHIKLMNPQRSTVWY-'


def _get_standard_operon():
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
    return operon


def test_calculate_adjusted_operon_bounds():
    operon = _get_standard_operon()
    offset, length = calculate_adjusted_operon_bounds(operon)
    assert offset == 12
    assert length == 1188


def test_count_results():
    csv_text = [
            'fail exclude:cas3 max-distance-to-anything:transposase-500 min-distance-to-anything:transposase-1',
            'fail require:transposase max-distance-to-anything:transposase-500 min-distance-to-anything:transposase-1',
            'pass',
            'fail exclude:cas3',
            'fail require:transposase',
            'fail exclude:cas3 max-distance-to-anything:transposase-500 min-distance-to-anything:transposase-1',
            'fail max-distance-to-anything:transposase-500',
            'fail exclude:cas3'
            ]
    unique_rule_violated, failed_rule_occurrences, rule_failure_counts = _count_results(csv_text)
    expected_urv = {'exclude:cas3': 2,
                    'require:transposase': 1,
                    'max-distance-to-anything:transposase-500': 1,
                    'min-distance-to-anything:transposase-1': 0}
    expected_fro = {'exclude:cas3': 4,
                    'max-distance-to-anything:transposase-500': 4,
                    'min-distance-to-anything:transposase-1': 3,
                    'require:transposase': 2}
    expected_rfc = {0: 1, 1: 4, 3: 3}
    assert unique_rule_violated == expected_urv
    assert failed_rule_occurrences == expected_fro
    assert rule_failure_counts == expected_rfc


@composite
def random_feature(draw):
    name = draw(text(name_characters, min_size=3, max_size=8))
    start = draw(integers(min_value=0, max_value=100000))
    end = draw(integers(min_value=0, max_value=100000))
    sequence = draw(text(sequence_characters, min_size=50))
    accession = draw(text(name_characters, min_size=8, max_size=16))
    e_val = draw(floats(allow_nan=False, allow_infinity=False))
    strand = draw(sampled_from([-1, 1, None]))
    return Feature(name, (start, end), f'lcl|{start}|{end}|1|{strand}', strand, accession, e_val, '', sequence)


@pytest.mark.slow
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_valid_operon(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', 0, 100000, ls)
    rs = RuleSet().same_orientation()
    rs.evaluate(operon)
    assert True


def test_not_same_orientation():
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (600, 410), 'lcl|410|600|1|-1', -1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
    rs = RuleSet().same_orientation()
    result = rs.evaluate(operon)
    assert not result.is_passing


def test_same_orientation():
    operon = _get_standard_operon()
    rs = RuleSet().same_orientation()
    result = rs.evaluate(operon)
    assert result.is_passing


@pytest.mark.parametrize("feature_names,expected", [
    (['cas1', 'cas2', 'cas4'], True),
    (['cas2', 'cas4'], True),
    (['cas1', 'cas2', 'cas4', 'cas5'], False),
    ])
def test_contains_features(feature_names, expected):
    operon = _get_standard_operon()
    assert _contains_features(operon, feature_names) is expected


def test_contains_any_set_of_serialization():
    rs = RuleSet().contains_any_set_of_features([['cas1', 'cas2'],
                                                 ['cas2', 'cas3', 'cas4']])
    rule = rs._rules[0]
    expected = "contains-any-set-of-features:cas1-cas2|cas2-cas3-cas4"
    assert str(rule) == expected


@pytest.mark.parametrize("sets,expected", [
    ([['cas1', 'cas2', 'cas4'], ['cas3', 'cas5']], True),
    ([['cas1', 'cas2', 'cas4'], ['cas1', 'cas2']], True),
    ([['cas1', 'cas2', 'cas3', 'cas4'], ['cas3', 'cas5']], False),
    ([['cas1', 'cas2', 'cas3', 'cas4']], False),
    ])
def test_contains_any_set_of_features(sets, expected):
    operon = _get_standard_operon()
    rs = RuleSet().contains_any_set_of_features(sets)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


@pytest.mark.parametrize("f1,f2,expected", [
    ('cas1', 'cas700', True),
    ('cas2', 'cas700', True),
    ('cas700', 'cas1', True),
    ('cas700', 'cas2', True),
    ('cas1', 'cas2', False),
    ('cas2', 'cas1', False),
    ('cas700', 'cas800', False),
    ])
def test_contains_exactly_one_of(f1, f2, expected):
    operon = _get_standard_operon()
    rs = RuleSet().contains_exactly_one_of(f1, f2)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


@pytest.mark.parametrize('distance,expected', [
    (0, True),
    (50, True),
    (99, True),
    (100, False),
    (1000, False),
    (10000, False)
    ])
def test_at_least_n_bp_from_anything(distance: int, expected: bool):
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('transposase', (700, 800), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('cas4', (920, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
    rs = RuleSet().at_least_n_bp_from_anything('transposase', distance)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


@pytest.mark.parametrize('distance,expected', [
    (0, False),
    (50, False),
    (99, False),
    (100, True),
    (1000, True),
    (10000, True)
    ])
def test_at_most_n_bp_from_anything(distance: int, expected: bool):
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('transposase', (700, 800), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('cas4', (920, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
    rs = RuleSet().at_most_n_bp_from_anything('transposase', distance)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


def test_serialize_results_fail():
    operon = _get_standard_operon()
    rs = RuleSet() \
        .exclude('cas3') \
        .require('cas12a')
    result = rs.evaluate(operon)
    actual = "\n".join(_serialize_results(rs, [result]))
    expected = "# exclude:cas3,require:cas12a\nQCDRTU,0..3400,fail require:cas12a"
    assert actual == expected


@pytest.mark.parametrize('gene1_start,gene1_end,gene2_start,gene2_end', [
    (12, 400, 410, 600),
    (410, 600, 12, 400),
    (400, 12, 410, 600),
    (12, 400, 600, 410),
    (410, 600, 400, 12),
    ])
def test_feature_distance(gene1_start, gene1_end, gene2_start, gene2_end):
    f1 = Feature('cas1', (gene1_start, gene1_end), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER')
    f2 = Feature('cas2', (gene2_start, gene2_end), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')
    distance = _feature_distance(f1, f2)
    distance_reverse = _feature_distance(f2, f1)
    assert distance == 10
    assert distance == distance_reverse


def test_feature_distance_overlapping():
    f1 = Feature('cas1', (100, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER')
    f2 = Feature('cas2', (200, 300), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')
    distance = _feature_distance(f1, f2)
    distance_reverse = _feature_distance(f2, f1)
    assert distance == 0
    assert distance == distance_reverse


@pytest.mark.parametrize('f1name,f2name,expected', [
    ('cas1', 'cas2', True),
    ('cas2', 'cas1', True),
    ('cas1', 'cas77', False),
    ('cas77', 'cas2', False),
    ('cas77', 'cas88', False),
    ('cas77', 'cas1', False),
    ('cas77', 'cas2', False)])
def test_max_distance_missing(f1name, f2name, expected):
    genes = [
            Feature('cas1', (0, 300), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (310, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            ]
    operon = Operon('contig', 0, 1000, genes)
    assert _max_distance(operon, f1name, f2name, 100) is expected


@pytest.mark.parametrize('gene1_start,gene1_end,gene2_start,gene2_end,distance_bp,expected', [
    (12, 400, 410, 600, 20, True),
    (410, 600, 12, 400, 20, True),
    (400, 12, 410, 600, 20, True),
    (12, 400, 600, 410, 20, True),
    (410, 600, 400, 12, 20, True),
    (12, 400, 410, 600, 5, False),
    (410, 600, 12, 400, 5, False),
    (400, 12, 410, 600, 5, False),
    (12, 400, 600, 410, 5, False),
    (410, 600, 400, 12, 5, False),
    ])
def test_max_distance(gene1_start, gene1_end, gene2_start, gene2_end, distance_bp, expected):
    genes = [
            Feature('cas1', (gene1_start, gene1_end), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (gene2_start, gene2_end), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            ]
    operon = Operon('contig', 0, 1000, genes)
    rs = RuleSet().max_distance('cas1', 'cas2', distance_bp)
    result = rs.evaluate(operon)
    assert result.is_passing is expected
