from operon_analyzer.genes import Feature, Operon
from operon_analyzer.visualize import calculate_adjusted_operon_bounds, create_operon_figure
from operon_analyzer.overview import _count_results
from operon_analyzer.parse import _parse_feature
from operon_analyzer.rules import RuleSet, FilterSet
import pytest
from matplotlib.text import Text


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


def _get_standard_operon():
    genes = [
            Feature('cas1', (12, 400), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', 0, 3400, genes)
    return operon


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


def _find_plotted_features(ax):
    # determines the names of all the features plotted in a dna_features_viewer plot
    features = set()
    for child in ax.properties()['children']:
        # we check if child._text is empty since there are blank text boxes for some reason
        if type(child) == Text and child._text:
            features.add(child._text)
    return features


def test_calculate_adjusted_operon_bounds():
    operon = _get_standard_operon()
    offset, length = calculate_adjusted_operon_bounds(operon)
    assert offset == 12
    assert length == 1188


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
