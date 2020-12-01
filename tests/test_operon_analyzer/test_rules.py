import random
import itertools
import pytest
import string
from hypothesis.strategies import composite, text, integers, sampled_from, floats, lists
from hypothesis import given, settings
from operon_analyzer import rules
from operon_analyzer.genes import Feature, Operon
from typing import List


def _get_repositionable_operon(s1, e1, s2, e2, s3, e3, s4, e4, arraystart, arrayend):
    genes = [
            Feature('cas1', (min(s1, e1), max(s1, e1)), 'lcl|12|400|1|-1', 1 if e1 > s1 else -1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (min(s2, e2), max(s2, e2)), 'lcl|410|600|1|-1', 1 if e2 > s2 else -1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('transposase', (min(s3, e3), max(s3, e3)), 'lcl|620|1200|1|-1', 1 if e3 > s3 else -1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('tnsA', (min(s4, e4), max(s4, e4)), 'lcl|620|1200|1|-1', 1 if e4 > s4 else -1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MTNSA'),
            Feature('CRISPR array', (arraystart, arrayend), '', None, '', None, 'CRISPR array with some repeats', 'ACGTTGATATTTATAGCGCA'),
            ]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, max(s1, s2, s3, s4, arraystart, e1, e2, e3, e4, arrayend), genes)
    return operon


def _get_standard_operon():
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    return operon


def _get_group_operon():
    genes = [
            Feature('gwrN', (0, 99), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('transposase', (100, 200), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas1', (210, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('mrtK', (1300, 1500), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    return Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)


def test_rule_group():
    operon = _get_group_operon()
    for feature_names in itertools.permutations(['transposase', 'cas1', 'cas2', 'cas4']):
        assert rules._contains_group(operon, feature_names, 20, False)


def test_rule_group_case_insensitive():
    operon = _get_group_operon()
    for feature_names in itertools.permutations(['Transposase', 'Cas1', 'Cas2', 'CAS4']):
        assert rules._contains_group(operon, feature_names, 20, False)


def test_rule_group_case_insensitive2():
    operon = _get_group_operon()
    for feature_names in itertools.permutations(['Cas2', 'CAS4', 'mrtk']):
        assert rules._contains_group(operon, feature_names, 120, False)


def test_rule_group_missing():
    operon = _get_group_operon()
    for feature_names in itertools.permutations(['transposase', 'cas9000', 'cas2', 'cas4']):
        assert not rules._contains_group(operon, feature_names, 20, False)


def test_rule_minimum_size():
    operon = _get_standard_operon()
    assert rules._minimum_size(operon, 'cas1', 350, False, False)
    assert not rules._minimum_size(operon, 'cas1', 400, False, False)


def test_rule_minimum_size_regex():
    operon = _get_standard_operon()
    assert rules._minimum_size(operon, r'cas\d', 100, False, True)
    assert rules._minimum_size(operon, r'cas\d', 100, True, True)
    assert rules._minimum_size(operon, r'cas\d', 300, False, True)
    assert not rules._minimum_size(operon, r'cas\d', 300, True, True)
    assert not rules._minimum_size(operon, r'cas\d', 700, False, True)


def test_rule_maximum_size():
    operon = _get_standard_operon()
    assert rules._maximum_size(operon, 'cas1', 450, False, False)
    assert not rules._maximum_size(operon, 'cas1', 380, False, False)


def test_rule_maximum_size_regex():
    operon = _get_standard_operon()
    assert rules._maximum_size(operon, r'cas\d', 200, False, True)
    assert rules._maximum_size(operon, r'cas\d', 700, True, True)
    assert rules._maximum_size(operon, r'cas\d', 400, False, True)
    assert not rules._maximum_size(operon, r'cas\d', 400, True, True)
    assert not rules._maximum_size(operon, r'cas\d', 100, False, True)


def test_rule_group_distance():
    operon = _get_group_operon()
    for feature_names in itertools.permutations(['transposase', 'cas1', 'cas2', 'cas4']):
        assert not rules._contains_group(operon, feature_names, 5, False)


def test_rule_group_same_orientation():
    operon = _get_group_operon()
    # all features have the same orientation in this operon to begin with
    for feature_names in itertools.permutations(['transposase', 'cas1', 'cas2', 'cas4']):
        assert rules._contains_group(operon, feature_names, 20, True)

    # ensure strands of unrelated Features aren't tested for orientation
    operon._features[0].strand = -1
    for feature_names in itertools.permutations(['transposase', 'cas1', 'cas2', 'cas4']):
        assert rules._contains_group(operon, feature_names, 20, True)

    # revert the strand of the unrelated Feature and then flip one that does matter
    operon._features[0].strand = 1
    operon._features[3].strand = -1
    for feature_names in itertools.permutations(['transposase', 'cas1', 'cas2', 'cas4']):
        assert not rules._contains_group(operon, feature_names, 20, True)


def _not_a_good_group_of_cas1_cas2_and_cas3_genes():
    genes = [
            Feature('cas1', (0, 70), 'lcl|0|70|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (100, 200), 'lcl|100|200|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas3', (240, 280), 'lcl|280|240|1|-1', -1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            ]
    return genes


def _a_good_group_of_cas1_cas2_and_cas3_genes():
    genes = [
            Feature('cas1', (1410, 1610), 'lcl|1410|1610|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas2', (1620, 2200), 'lcl|1620|2200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('cas3', (2210, 2500), 'lcl|2210|2500|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE')
            ]
    return genes


def test_contains_group_premature_return():
    # tests an edge case where the operon has two groups that contain all required features,
    # but only the second group has a passing combination of features, distances, and orientations
    # Previously, the first (non-valid) group would have caused the rule to return False prematurely

    # first, assert that the non-valid group really doesn't pass, even though it has the Features we want
    bad_operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, _not_a_good_group_of_cas1_cas2_and_cas3_genes())
    for feature_names in itertools.permutations(['cas1', 'cas2', 'cas3']):
        assert not rules._contains_group(bad_operon, feature_names, 20, False)
        assert not rules._contains_group(bad_operon, feature_names, 50, True)
    
    # now, append a "passing" group of features. The same rules should pass now
    genes = _not_a_good_group_of_cas1_cas2_and_cas3_genes() + \
            [Feature('cas4', (1200, 1300), 'lcl|1200|1300|1|-1', -1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER')] + \
            _a_good_group_of_cas1_cas2_and_cas3_genes()
    good_operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    for feature_names in itertools.permutations(['cas1', 'cas2', 'cas3']):
        assert rules._contains_group(good_operon, feature_names, 20, False)
        assert rules._contains_group(good_operon, feature_names, 50, True)


def test_case_insensitivity():
    genes = [
            Feature('Cas1', (12, 400), 'lcl|12|400|1|1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('Cas2', (410, 600), 'lcl|410|600|1|-1', -1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('Cas4', (620, 1200), 'lcl|620|1200|1|1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    assert rules._contains_features(operon, ['cas1', 'cas2'])
    assert rules._at_least_n_bp_from_anything(operon, 'cas1', 5)
    assert rules._at_least_n_bp_from_anything(operon, 'cas1', 5, regex=True)
    assert rules._at_most_n_bp_from_anything(operon, 'cas1', 100)
    assert rules._at_most_n_bp_from_anything(operon, 'cas1', 100, regex=True)
    assert rules._contains_any_set_of_features(operon, [['cas1', 'cas2']])
    assert rules._contains_at_least_n_features(operon, ['cas1', 'cas2', 'cas4'], 1, False)
    assert rules._contains_exactly_one_of(operon, 'cas1', 'cas9')
    assert rules._contains_features(operon, ['cas1', 'cas2', 'cas4'])
    assert rules._exclude(operon, 'cas9', False)
    assert rules._exclude(operon, 'cas9', True)
    assert rules._max_distance(operon, 'cas1', 'cas2', 20, True, False)
    assert not rules._max_distance(operon, 'cas1', 'cas2', 5, True, False)
    assert rules._require(operon, 'cas1', False)
    assert rules._require(operon, 'cas1', True)
    assert not rules._require(operon, 'cas9', False)
    assert not rules._require(operon, 'cas9', True)
    assert rules._same_orientation(operon, exceptions=['cas2'])
    assert not rules._same_orientation(operon, exceptions=[])


def test_ignore_filtered_features_in_rule():
    op = _get_standard_operon()
    rs = rules.RuleSet().exclude('cas4')
    assert not rs.evaluate(op).is_passing
    fs = rules.FilterSet().must_be_within_n_bp_of_anything(15)
    fs.evaluate(op)
    assert rs.evaluate(op).is_passing


def test_exclude_regex():
    op = _get_standard_operon()
    rs = rules.RuleSet().exclude(r'cas\d', True)
    assert not rs.evaluate(op).is_passing


def test_require_regex():
    op = _get_standard_operon()
    rs = rules.RuleSet().require(r'cas\d', True)
    assert rs.evaluate(op).is_passing


@pytest.mark.parametrize('f1,f2,distance,expected', [
    ('cas2', 'cas', 10, True),
    ('cas4', 'cas', 19, False),
    ('cas4', 'cas', 20, True),
    ('cas', 'cas', 10, True),
    ('cas', 'cas', 9, False),
    ('cas2', 'cas4', 19, False),
    ('cas2', 'cas4', 20, True),
    (r'cas\d+', r'cas\d+', 10, True),
    (r'cas\d\d', r'cas\d\d', 10, False),
    (r'.*?\d', r'cas\d', 10, True),
    ])
def test_max_distance_closest_regex(f1, f2, distance, expected):
    op = _get_standard_operon()
    rs = rules.RuleSet().max_distance(f1, f2, distance, closest_pair_only=True, regex=True)
    result = rs.evaluate(op)
    assert result.is_passing is expected


def test_at_least_n_bp_from_anything_regex():
    op = _get_standard_operon()
    rs = rules.RuleSet().at_least_n_bp_from_anything(r'cas\d', 10, regex=True)
    assert rs.evaluate(op).is_passing


@pytest.mark.parametrize('name,distance,expected', [
    (r'cas\d', 20, True),
    (r'cas\d', 19, False),
    (r'cas', 20, True),
    (r'cas', 19, False),
    (r'lulz', 10, False),
    ])
def test_at_most_n_bp_from_anything_regex(name, distance, expected):
    op = _get_standard_operon()
    rs = rules.RuleSet().at_most_n_bp_from_anything(name, distance, regex=True)
    assert rs.evaluate(op).is_passing is expected


@pytest.mark.parametrize('f1,f2,expected', [
    (r'cas', r'lulz', True),
    (r'lulz', r'cas', True),
    (r'histone', r'lulz', False),
    (r'cas', r'cas\d', False),
    ])
def test_contains_exactly_one_of_regex(f1, f2, expected):
    op = _get_standard_operon()
    rs = rules.RuleSet().contains_exactly_one_of(f1, f2, regex=True)
    assert rs.evaluate(op).is_passing is expected


@pytest.mark.parametrize('feature,expected_count', [
    ('cas', 3),
    ('cas4', 1),
    ('lulz', 0),
    ])
def test_regex_rule(feature, expected_count):
    op = _get_standard_operon()
    features = op.get(feature, regex=True)
    assert len(features) == expected_count


def test_filterset_within_n_bp_anything():
    operon = _get_standard_operon()
    fs = rules.FilterSet().must_be_within_n_bp_of_anything(10)
    fs.evaluate(operon)
    names = list(operon.feature_names)
    assert 'cas4' not in names
    assert 'cas1' in names
    assert 'cas2' in names


def test_filterset_within_n_bp_anything_one_feature():
    genes = [Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER')]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    fs = rules.FilterSet().must_be_within_n_bp_of_anything(10)
    fs.evaluate(operon)
    names = list(operon.feature_names)
    assert names == ['cas1']


def test_filterset_within_n_bp_of_feature():
    operon = _get_standard_operon()
    fs = rules.FilterSet().must_be_within_n_bp_of_feature('cas2', 10)
    fs.evaluate(operon)
    names = list(operon.feature_names)
    assert 'cas4' not in names
    assert 'cas1' in names
    assert 'cas2' in names


def test_custom_rule():
    operon = _get_standard_operon()
    rule = rules.Rule('require', rules._require, 'cas1', None)
    rs = rules.RuleSet().custom(rule)
    result = rs.evaluate(operon)
    assert result.is_passing

    rule = rules.Rule('require', rules._require, 'cas88', None)
    rs = rules.RuleSet().custom(rule)
    result = rs.evaluate(operon)
    assert not result.is_passing


def test_filterset_within_n_bp_of_feature_only_one_feature():
    genes = [Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER')]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    fs = rules.FilterSet().must_be_within_n_bp_of_feature('cas1', 10)
    fs.evaluate(operon)
    names = list(operon.feature_names)
    assert 'cas1' in names


def test_at_most_n_bp_single_feature():
    """ Ensure we don't crash when a single Feature is present. """
    genes = [Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER')]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    rs = rules.RuleSet().at_most_n_bp_from_anything('cas1', 50)
    rs.evaluate(operon)
    assert True


def test_filter_overlap_reason_text():
    bit_scores = [100.0, 100.0, 100.0, 200.0]
    positions = [0, 100, 101, 200, 201, 300, 210, 300, 400, 500]
    operon = _get_repositionable_operon(*positions)
    for gene, bit_score in zip(operon.all_genes, bit_scores):
        gene.bit_score = bit_score
    fs = rules.FilterSet().pick_overlapping_features_by_bit_score(0.8)
    fs.evaluate(operon)
    for feature in operon.all_features:
        if feature.name == 'transposase':
            assert feature.ignored_reasons == ['overlaps-tnsA:0.8']
            break
    else:
        # make sure we ran the above assertion and broke out of the loop
        assert False


@pytest.mark.parametrize('positions,bit_scores,expected,threshold', [
    ([0, 100, 101, 200, 201, 300, 210, 300, 400, 500], [100.0, 100.0, 100.0, 200.0], [False, False, True, False], 0.8),
    ([0, 100, 101, 200, 201, 300, 210, 300, 400, 500], [100.0, 100.0, 200.0, 100.0], [False, False, False, True], 0.8),
    ([0, 100, 101, 200, 201, 300, 310, 400, 500, 600], [100.0, 100.0, 200.0, 100.0], [False, False, False, False], 0.8),
    ([0, 100, 101, 200, 201, 300, 211, 300, 400, 500], [100.0, 100.0, 100.0, 200.0], [False, False, False, False], 0.99),
    ([0, 100, 101, 200, 201, 300, 211, 300, 400, 500], [100.0, 100.0, 100.0, 100.0], [False, False, False, False], 0.8),
    ([0, 100, 101, 200, 201, 300, 201, 300, 400, 500], [100.0, 100.0, 105.0, 100.0], [False, False, False, True], 0.8),
    ([1, 200, 141, 180, 141, 300, 500, 600, 700, 800], [444.0, 100.0, 444.0, 100.0], [False, True, False, False], 0.8),
    ])
def test_pick_overlapping_features_by_bit_score(positions: List[int],
                                                bit_scores: List[float],
                                                expected: List[bool],
                                                threshold: float,
                                                ):
    operon = _get_repositionable_operon(*positions)
    for gene, bit_score in zip(operon.all_genes, bit_scores):
        gene.bit_score = bit_score
    fs = rules.FilterSet().pick_overlapping_features_by_bit_score(threshold)
    fs.evaluate(operon)
    actual = [bool(feature.ignored_reasons) for feature in operon.all_genes]
    assert expected == actual


def _get_standard_operon_with_overlapping_feature():
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER', 152),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR', 143),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE', 546),
            Feature('cas5', (630, 1211), 'lcl|630|1211|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'an overlapping gene', 'MLAWPVTLE', 93.7),
            ]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    return operon


def test_pick_overlapping_features_by_bit_score_2():
    expected = [False, False, False, True]
    operon = _get_standard_operon_with_overlapping_feature()
    fs = rules.FilterSet().pick_overlapping_features_by_bit_score(0.8)
    fs.evaluate(operon)
    actual = [bool(feature.ignored_reasons) for feature in operon.all_genes]
    assert expected == actual


@pytest.mark.parametrize('fstart,fend,ostart,oend,expected', [
    (1, 100, 51, 150, 0.5),
    (51, 150, 1, 100, 0.5),
    (0, 100, 200, 300, None),
    (200, 300, 0, 100, None),
    (1, 2, 2, 3, 0.5),
    (2, 3, 1, 2, 0.5),
    (0, 5, 3, 8, 0.5),
    (3, 8, 0, 5, 0.5),
    (1, 100, 91, 200, 0.1),
    (1, 100, 11, 200, 0.9),
    (1, 100, 101, 200, None),
    (101, 200, 1, 100, None),
    ])
def test_calculate_overlap(fstart, fend, ostart, oend, expected):
    feature = Feature("tnsA", (fstart, fend), "", 1, "", 0.001, "", "MRTK")
    other_feature = Feature("transposase", (ostart, oend), "", 1, "", 0.001, "", "MGWRN")
    overlap = rules._calculate_overlap(feature, other_feature)
    if expected is None:
        assert overlap is None
    else:
        assert pytest.approx(overlap) == expected


@pytest.mark.parametrize('feature_name,expected', [
    ('cas2', True),
    ('cas1', False)
    ])
def test_multicopy_feature(feature_name, expected):
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('cas2', (1220, 1300), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas1', (2000, 2200), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            ]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    rs = rules.RuleSet().at_most_n_bp_from_anything(feature_name, 25)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


name_characters = string.ascii_lowercase + string.ascii_uppercase + string.digits

sequence_characters = 'ACDEFGHIKLMNPQRSTVWY-'

@composite
def random_feature(draw):
    name = draw(text(name_characters, min_size=3, max_size=8))
    start = draw(integers(min_value=0, max_value=99990))
    end = draw(integers(min_value=start+3, max_value=100000))
    sequence = draw(text(sequence_characters, min_size=50))
    accession = draw(text(name_characters, min_size=8, max_size=16))
    e_val = draw(floats(allow_nan=False, allow_infinity=False))
    strand = draw(sampled_from([-1, 1, None]))
    return Feature(name, (start, end), f'lcl|{start}|{end}|1|{strand}', strand, accession, e_val, '', sequence)


@pytest.mark.proptest
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_all_fixed_rules(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', '/tmp/dna.fasta', 0, 100000, ls)
    rs = rules.RuleSet().same_orientation() \
                  .at_least_n_bp_from_anything('cas1', 50) \
                  .at_most_n_bp_from_anything('cas2', 100) \
                  .contains_any_set_of_features([['cas5', 'cas6', 'cas7'],
                                                 ['cas12a'],
                                                 ['cas4', 'cas10']]) \
                  .contains_exactly_one_of('cas1', 'cas13') \
                  .exclude('lulz') \
                  .max_distance('f1', 'f2', 100) \
                  .require('woo') \
                  .same_orientation()

    rs.evaluate(operon)
    assert True


@pytest.mark.proptest
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_at_least_n_bp_from_anything(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', '/tmp/dna.fasta', 0, 100000, ls)
    random_feature = random.choice(list(operon.feature_names))
    rs = rules.RuleSet().at_least_n_bp_from_anything(random_feature, 50)
    rs.evaluate(operon)
    assert True


@pytest.mark.proptest
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_at_most_n_bp_from_anything(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', '/tmp/dna.fasta', 0, 100000, ls)
    random_feature = random.choice(list(operon.feature_names))
    rs = rules.RuleSet().at_most_n_bp_from_anything(random_feature, 50)
    rs.evaluate(operon)
    assert True


@pytest.mark.proptest
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_contains_any_set(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', '/tmp/dna.fasta', 0, 100000, ls)
    random_features1 = random.sample(list(operon.feature_names), random.randint(1, len(operon)))
    random_features2 = random.sample(list(operon.feature_names), random.randint(1, len(operon)))
    rs = rules.RuleSet().contains_any_set_of_features([random_features1, random_features2])
    rs.evaluate(operon)
    assert True


@pytest.mark.proptest
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_contains_exactly_one_of(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', '/tmp/dna.fasta', 0, 100000, ls)
    random_feature = random.choice(list(operon.feature_names))
    random_feature2 = random.choice([random.choice(list(operon.feature_names)), "lulz"])
    rs = rules.RuleSet().contains_exactly_one_of(random_feature, random_feature2)
    rs.evaluate(operon)
    assert True


@pytest.mark.proptest
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_exclude(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', '/tmp/dna.fasta', 0, 100000, ls)
    random_feature = random.choice([random.choice(list(operon.feature_names)), "lulz"])
    rs = rules.RuleSet().exclude(random_feature)
    rs.evaluate(operon)
    assert True


@pytest.mark.proptest
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_max_distance(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', '/tmp/dna.fasta', 0, 100000, ls)
    random_feature = random.choice(list(operon.feature_names))
    random_feature2 = random.choice(list(operon.feature_names))
    distance_bp = random.randint(0, 99999)
    rs = rules.RuleSet().max_distance(random_feature, random_feature2, distance_bp)
    rs.evaluate(operon)
    assert True


@pytest.mark.proptest
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_require(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', '/tmp/dna.fasta', 0, 100000, ls)
    random_feature = random.choice(list(operon.feature_names))
    rs = rules.RuleSet().require(random_feature)
    rs.evaluate(operon)
    assert True


@pytest.mark.proptest
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_same_orientation(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', '/tmp/dna.fasta', 0, 100000, ls)
    count = random.randint(1, len(operon))
    exceptions = random.choice([random.sample(list(operon.feature_names), count), None])
    rs = rules.RuleSet().same_orientation(exceptions=exceptions)
    rs.evaluate(operon)
    assert True


def test_not_same_orientation():
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', -1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    rs = rules.RuleSet().same_orientation()
    result = rs.evaluate(operon)
    assert not result.is_passing


def test_same_orientation():
    operon = _get_standard_operon()
    rs = rules.RuleSet().same_orientation()
    result = rs.evaluate(operon)
    assert result.is_passing


@pytest.mark.parametrize('exceptions,expected', [
        (['cas2', 'CRISPR array'], True),
        (['cas2'], False),
        (['CRISPR array'], False),
        ([], False),
        (None, False)
    ])
def test_same_orientation_with_exceptions(exceptions, expected):
    positions = [0, 100, 701, 600, 201, 300, 210, 300, 400, 500]
    operon = _get_repositionable_operon(*positions)
    rs = rules.RuleSet().same_orientation(exceptions=exceptions)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


@pytest.mark.parametrize("feature_names,expected", [
    (['cas1', 'cas2', 'cas4'], True),
    (['cas2', 'cas4'], True),
    (['cas1', 'cas2', 'cas4', 'cas5'], False),
    ])
def test_contains_features(feature_names, expected):
    operon = _get_standard_operon()
    assert rules._contains_features(operon, feature_names) is expected


def test_contains_any_set_of_serialization():
    rs = rules.RuleSet().contains_any_set_of_features([['cas1', 'cas2'],
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
    rs = rules.RuleSet().contains_any_set_of_features(sets)
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
    rs = rules.RuleSet().contains_exactly_one_of(f1, f2)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


@pytest.mark.parametrize('distance,expected', [
    (0, True),
    (50, True),
    (99, True),
    (100, True),
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
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    rs = rules.RuleSet().at_least_n_bp_from_anything('transposase', distance)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


def test_at_least_n_bp_from_anything_no_distances():
    genes = [Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER')]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    rs = rules.RuleSet().at_least_n_bp_from_anything('cas1', 5)
    result = rs.evaluate(operon)
    assert result.is_passing is True


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
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    rs = rules.RuleSet().at_most_n_bp_from_anything('transposase', distance)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


@pytest.mark.parametrize('gene1_start,gene1_end,gene2_start,gene2_end', [
    (12, 400, 410, 600),
    (410, 600, 12, 400),
    ])
def test_feature_distance(gene1_start, gene1_end, gene2_start, gene2_end):
    f1 = Feature('cas1', (gene1_start, gene1_end), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER')
    f2 = Feature('cas2', (gene2_start, gene2_end), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')
    distance = rules._feature_distance(f1, f2)
    distance_reverse = rules._feature_distance(f2, f1)
    assert distance == 10
    assert distance == distance_reverse


def test_feature_distance_overlapping():
    f1 = Feature('cas1', (100, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER')
    f2 = Feature('cas2', (200, 300), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')
    distance = rules._feature_distance(f1, f2)
    distance_reverse = rules._feature_distance(f2, f1)
    assert distance == 0
    assert distance == distance_reverse


@pytest.mark.parametrize('f1name,f2name,expected', [
    ('cas1', 'cas1', False),
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
    operon = Operon('contig', '/tmp/dna.fasta', 0, 1000, genes)
    rs = rules.RuleSet().max_distance(f1name, f2name, 100)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


@pytest.mark.parametrize('gene1_start,gene1_end,gene2_start,gene2_end,distance_bp,expected', [
    (12, 400, 410, 600, 20, True),
    (410, 600, 12, 400, 20, True),
    (12, 400, 410, 600, 5, False),
    (410, 600, 12, 400, 5, False),
    ])
def test_max_distance(gene1_start, gene1_end, gene2_start, gene2_end, distance_bp, expected):
    genes = [
            Feature('Cas1', (gene1_start, gene1_end), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('Cas2', (gene2_start, gene2_end), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            ]
    operon = Operon('contig', '/tmp/dna.fasta', 0, 1000, genes)
    rs = rules.RuleSet().max_distance('cas1', 'cas2', distance_bp)
    result = rs.evaluate(operon)
    assert result.is_passing is expected
