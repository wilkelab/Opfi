from operon_analyzer.genes import Feature, Operon
from operon_analyzer.rules import RuleSet, _feature_distance, _max_distance, _contains_features
from operon_analyzer.analyze import _serialize_results
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


@pytest.mark.parametrize('distance,expected', [
    (0, True),
    (50, True),
    (99, True),
    (100, False),
    (1000, False),
    (10000, False)
    ])
def test_min_distance_to_anything(distance: int, expected: bool):
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('transposase', (700, 800), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('cas4', (920, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
    rs = RuleSet().min_distance_to_anything('transposase', distance)
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
def test_max_distance_to_anything(distance: int, expected: bool):
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('transposase', (700, 800), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            Feature('cas4', (920, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
    rs = RuleSet().max_distance_to_anything('transposase', distance)
    result = rs.evaluate(operon)
    assert result.is_passing is expected


def test_serialize_results_fail():
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
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


# def test_parse_input_crispr_repeat():
#     data = ['STDS9ITBSE', '1009983..1031284', 'CRISPR array', '1017560..1017914', '', '', '', '', 'Copies: 6, Repeat: 30, Spacer: 34', '-------------------GATAAACATTAACATAGGATGTATTGAAAC-']


# def test_parse_input_cas_gene():
# data = ['CRISPR-transposases', '1396535..1417167', 'cas11', '1412833..1412048', 'lcl|1412833|1412048|2|-1', '-1', 'UniRef50_UPI0009FEB4C2', '0', 'type I-E CRISPR-associated protein Cse2/CasB n=1 Tax=Actinobaculum sp. oral taxon 183 TaxID=712888 RepID=UPI0009FEB4C2', 'MVRHRPKQSPSYIYHFPTSERERSIVTTVNEIIKDEKPLRKRKRRNLSPIGQKIDCKISCLQKGYLSEDSRKQARARADLANLRRGLTAGPGERVEIWHLTQVDVSDNAPDEPTREEFAVHVSMTLYAAHQQSRTKPMHRPAEGLGHAAHSVVGYGDDENPSARARFDALVMSSTPRELRRHLRSFVSLLRAKEIPLDYGMLVDDIVCFQRPGGAKAVRRHWSRQYYDFSSTDGESSETDSTAEDICSENSLHNSLRNTKE']
