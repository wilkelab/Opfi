from operon_analyzer.genes import Feature, Operon
from operon_analyzer.rules import RuleSet, _feature_distance, _max_distance, _contains_features, FilterSet, _calculate_overlap, _pick_overlapping_features_by_bit_score
from operon_analyzer.visualize import calculate_adjusted_operon_bounds, create_operon_figure
from operon_analyzer.overview import _count_results
from operon_analyzer.parse import _parse_feature
import pytest
from hypothesis.strategies import composite, text, integers, sampled_from, floats, lists
from hypothesis import given, settings
from typing import List
import string
from matplotlib.text import Text
import random


name_characters = string.ascii_lowercase + string.ascii_uppercase + string.digits

sequence_characters = 'ACDEFGHIKLMNPQRSTVWY-'


def test_at_most_n_bp_single_feature():
    """ Ensure we don't crash when a single Feature is present. """
    genes = [Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER')]
    operon = Operon('QCDRTU', 0, 3400, genes)
    rs = RuleSet().at_most_n_bp_from_anything('cas1', 50)
    rs.evaluate(operon)
    assert True


def _get_repositionable_operon(s1, e1, s2, e2, s3, e3, s4, e4, arraystart, arrayend):
    genes = [
            Feature('cas1', (s1, e1), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (s2, e2), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('transposase', (s3, e3), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('tnsA', (s4, e4), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MTNSA'),
            Feature('CRISPR array', (arraystart, arrayend), '', None, '', None, 'CRISPR array with some repeats', 'ACGTTGATATTTATAGCGCA'),
            ]
    operon = Operon('QCDRTU', 0, max(s1, s2, s3, s4, arraystart, e1, e2, e3, e4, arrayend), genes)
    return operon


def test_filter_overlap_reason_text():
    bit_scores = [100.0, 100.0, 100.0, 200.0]
    positions = [0, 100, 101, 200, 201, 300, 210, 300, 400, 500]
    operon = _get_repositionable_operon(*positions)
    for gene, bit_score in zip(operon.all_genes, bit_scores):
        gene.bit_score = bit_score
    fs = FilterSet().pick_overlapping_features_by_bit_score(0.8)
    fs.evaluate(operon)
    feature = operon.get_unique('transposase')
    assert feature.ignored_reasons == ['overlaps-tnsA:0.8']


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
    _pick_overlapping_features_by_bit_score(operon, 'overlaps-%s', threshold)
    actual = [bool(feature.ignored_reasons) for feature in operon.all_genes]
    assert expected == actual

def _get_standard_operon_with_overlapping_feature():
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER', 152),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR', 143),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE', 546),
            Feature('cas5', (630, 1211), 'lcl|630|1211|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'an overlapping gene', 'MLAWPVTLE', 93.7),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
    return operon

def test_pick_overlapping_features_by_bit_score_2():
    expected = [False, False, False, True]
    operon = _get_standard_operon_with_overlapping_feature()
    _pick_overlapping_features_by_bit_score(operon, 'overlaps-%s', 0.8)
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
    overlap = _calculate_overlap(feature, other_feature)
    if expected is None:
        assert overlap is None
    else:
        assert pytest.approx(overlap) == expected


@pytest.mark.parametrize('line,expected', [
    (('GDBD23958235',
      '327464..369995',
      'CRISPR array',
      '369885..369948',
      '',
      '',
      '',
      '',
      "Copies: 2, Repeat: 18, Spacer: 27",
      'AAGAAGGCTGCTAAGGTA'), "CRISPR array"),
    (('GBDB23958235',
     '534183..567749',
     'transposase',
     '545109..544183',
     'lcl|545109|544183|1|-1',
     '-1',
     'UniRef50_A0A1E3AYK0',
     '2.03206e-20',
     'nuclease family transposase n=26 Tax=Bacteria TaxID=2 RepID=A0A1E3AYK0_9FIRM',
     'EDKVFGMVMENKDFCKYLLEIIIPDLKIKKIDWLDKQVEINNSERK----NEAKEVRLDVLVTDHEGRVFNIEMQTTDQDDIGRRMRYYLSRLDLRYTLNKGKTYRNLKDAFIIFLCNFKPKKDDKFYESYHTYSDQDRSKQSQDGVTKIIINSQVSAEGQSEELKALAKLMNNEPVKLNKHFDYA-----QRRIKEINEDPEMREKIMLYETRMLEREQAAGKAGYEQ'), 'transposase'),
    ])
def test_parse_feature_name(line, expected):
    _, _, feature = _parse_feature(line)
    assert feature.name == expected



def _get_standard_operon():
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
    return operon


def _find_plotted_features(ax):
    # determines the names of all the features plotted in a dna_features_viewer plot
    features = set()
    for child in ax.properties()['children']:
        # we check if child._text is empty since there are blank text boxes for some reason
        if type(child) == Text and child._text:
            features.add(child._text)
    return features


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
    operon = Operon('QCDRTU', 0, 3400, genes)
    rs = RuleSet().at_most_n_bp_from_anything(feature_name, 25)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


def test_create_operon_figure_with_CRISPR_array():
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('CRISPR array', (1300, 1400), '', None, '', None, 'Copies: 2, Repeat: 61, Spacer: 59', 'CTCAAACTCTCACTCTGGCTCAATGAGTTAGACAAGCTCTCACTCTGACTCAAGGAATTAC'),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
    fs = FilterSet().must_be_within_n_bp_of_feature('cas2', 10000)
    rs = RuleSet().require('CRISPR array')
    fs.evaluate(operon)
    result = rs.evaluate(operon)
    assert result.is_passing
    ax = create_operon_figure(operon, True, None)
    features = _find_plotted_features(ax)
    assert features == set(['cas1', 'cas2', 'cas4', 'CRISPR array (2)'])


def test_create_operon_figure_with_ignored():
    operon = _get_standard_operon()
    fs = FilterSet().must_be_within_n_bp_of_feature('cas2', 10)
    fs.evaluate(operon)
    ax = create_operon_figure(operon, True, None)
    features = _find_plotted_features(ax)
    assert features == set(['cas1', 'cas2', 'cas4 (ignored)'])


def test_create_operon_figure_with_colors():
    operon = _get_standard_operon()
    fs = FilterSet().must_be_within_n_bp_of_feature('cas2', 10)
    fs.evaluate(operon)
    gene_colors = {'cas1': 'purple', 'cas2': 'green'}
    ax = create_operon_figure(operon, False, gene_colors)
    features = _find_plotted_features(ax)
    assert features == set(['cas1', 'cas2'])


def test_create_operon_figure():
    operon = _get_standard_operon()
    fs = FilterSet().must_be_within_n_bp_of_feature('cas2', 10)
    fs.evaluate(operon)
    ax = create_operon_figure(operon, False, None)
    features = _find_plotted_features(ax)
    assert features == set(['cas1', 'cas2'])


def test_filterset_within_n_bp_of_feature():
    operon = _get_standard_operon()
    fs = FilterSet().must_be_within_n_bp_of_feature('cas2', 10)
    fs.evaluate(operon)
    names = list(operon.feature_names)
    assert 'cas4' not in names
    assert 'cas1' in names
    assert 'cas2' in names


def test_filterset_within_n_bp_anything():
    operon = _get_standard_operon()
    fs = FilterSet().must_be_within_n_bp_of_anything(10)
    fs.evaluate(operon)
    names = list(operon.feature_names)
    assert 'cas4' not in names
    assert 'cas1' in names
    assert 'cas2' in names


def test_calculate_adjusted_operon_bounds():
    operon = _get_standard_operon()
    offset, length = calculate_adjusted_operon_bounds(operon)
    assert offset == 12
    assert length == 1188


def test_count_results():
    csv_text = [line.split(',') for line in [
            'fail,exclude:cas3,max-distance-to-anything:transposase-500,min-distance-to-anything:transposase-1',
            'fail,require:transposase,max-distance-to-anything:transposase-500,min-distance-to-anything:transposase-1',
            'pass',
            'fail,exclude:cas3',
            'fail,require:transposase',
            'fail,exclude:cas3,max-distance-to-anything:transposase-500,min-distance-to-anything:transposase-1',
            'fail,max-distance-to-anything:transposase-500',
            'fail,exclude:cas3']
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
    start = draw(integers(min_value=0, max_value=99990))
    end = draw(integers(min_value=start+3, max_value=100000))
    sequence = draw(text(sequence_characters, min_size=50))
    accession = draw(text(name_characters, min_size=8, max_size=16))
    e_val = draw(floats(allow_nan=False, allow_infinity=False))
    strand = draw(sampled_from([-1, 1, None]))
    return Feature(name, (start, end), f'lcl|{start}|{end}|1|{strand}', strand, accession, e_val, '', sequence)


@pytest.mark.slow
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_all_fixed_rules(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', 0, 100000, ls)
    rs = RuleSet().same_orientation() \
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


@pytest.mark.slow
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_at_least_n_bp_from_anything(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', 0, 100000, ls)
    random_feature = random.choice(list(operon.feature_names))
    rs = RuleSet().at_least_n_bp_from_anything(random_feature, 50)
    rs.evaluate(operon)
    assert True


@pytest.mark.slow
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_at_most_n_bp_from_anything(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', 0, 100000, ls)
    random_feature = random.choice(list(operon.feature_names))
    rs = RuleSet().at_most_n_bp_from_anything(random_feature, 50)
    rs.evaluate(operon)
    assert True


@pytest.mark.slow
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_contains_any_set(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', 0, 100000, ls)
    random_features1 = random.sample(list(operon.feature_names), random.randint(1, len(operon)))
    random_features2 = random.sample(list(operon.feature_names), random.randint(1, len(operon)))
    rs = RuleSet().contains_any_set_of_features([random_features1, random_features2])
    rs.evaluate(operon)
    assert True


@pytest.mark.slow
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_contains_exactly_one_of(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', 0, 100000, ls)
    random_feature = random.choice(list(operon.feature_names))
    random_feature2 = random.choice([random.choice(list(operon.feature_names)), "lulz"])
    rs = RuleSet().contains_exactly_one_of(random_feature, random_feature2)
    rs.evaluate(operon)
    assert True


@pytest.mark.slow
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_exclude(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', 0, 100000, ls)
    random_feature = random.choice([random.choice(list(operon.feature_names)), "lulz"])
    rs = RuleSet().exclude(random_feature)
    rs.evaluate(operon)
    assert True


@pytest.mark.slow
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_max_distance(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', 0, 100000, ls)
    random_feature = random.choice(list(operon.feature_names))
    random_feature2 = random.choice(list(operon.feature_names))
    distance_bp = random.randint(0, 99999)
    rs = RuleSet().max_distance(random_feature, random_feature2, distance_bp)
    rs.evaluate(operon)
    assert True


@pytest.mark.slow
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_require(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', 0, 100000, ls)
    random_feature = random.choice(list(operon.feature_names))
    rs = RuleSet().require(random_feature)
    rs.evaluate(operon)
    assert True


@pytest.mark.slow
@given(lists(random_feature(), min_size=1, max_size=20, unique=True))
@settings(max_examples=1000)
def test_random_operon_same_orientation(ls):
    """
    Just fuzzes evaluation with randomly-generated operons.
    All we're trying to do here is see if we can cause an uncaught exception.
    """
    operon = Operon('testoperon', 0, 100000, ls)
    count = random.randint(1, len(operon))
    exceptions = random.choice([random.sample(list(operon.feature_names), count), None])
    rs = RuleSet().same_orientation(exceptions=exceptions)
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


@pytest.mark.parametrize('feature_list,feature_count,must_be_unique,expected', [
    (['cas3', 'cas7', 'cas12'], 3, False, False),
    (['cas2', 'cas7', 'cas12'], 3, False, False),
    (['cas1', 'cas2', 'cas12'], 3, False, True),
    (['cas1', 'cas2', 'cas12'], 3, True, False),
    (['cas1', 'cas2', 'cas3', 'cas4', 'cas5'], 4, False, True),
    (['cas1', 'cas2', 'cas3', 'cas4', 'cas5'], 5, False, False),
    (['cas1', 'cas2', 'cas3', 'cas4', 'cas5'], 4, True, False)
    ])
def test_contains_at_least_n_features(feature_list, feature_count, must_be_unique, expected):
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas2', (650, 660), 'lcl|650|660|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAM'),
            Feature('transposase', (700, 800), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('cas4', (920, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
    rs = RuleSet().contains_at_least_n_features(feature_list, feature_count, must_be_unique)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


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
