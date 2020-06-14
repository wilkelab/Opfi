from matplotlib.text import Text
from operon_analyzer.genes import Feature, Operon
from operon_analyzer.rules import RuleSet, FilterSet
from operon_analyzer.visualize import calculate_adjusted_operon_bounds, \
                                      create_operon_figure, \
                                      build_image_filename, \
                                      build_operon_dictionary
from common import get_standard_operon
import pytest


@pytest.mark.parametrize('directory,expected', [
    (None, 'QCDRTU-0-3400.png'),
    ('images', 'images/QCDRTU-0-3400.png'),
    ('images/', 'images/QCDRTU-0-3400.png'),
    ])
def test_build_image_filename(directory: str, expected: str):
    assert build_image_filename(get_standard_operon(), directory) == expected


def _find_plotted_features(ax):
    # determines the names of all the features plotted in a dna_features_viewer plot
    features = set()
    for child in ax.properties()['children']:
        # we check if child._text is empty since there are blank text boxes for some reason
        if type(child) == Text and child._text:
            features.add(child._text)
    return features


def test_calculate_adjusted_operon_bounds_all_ignored():
    operon = get_standard_operon()
    for feature in operon:
        feature.ignore('')
    with pytest.raises(AssertionError):
        calculate_adjusted_operon_bounds(operon, False)
        assert True


def test_create_operon_figure_all_ignored():
    operon = get_standard_operon()
    for feature in operon:
        feature.ignore('')
    result = create_operon_figure(operon, False, None)
    assert result is None


def test_create_operon_figure_all_ignored_plot_ignored():
    operon = get_standard_operon()
    for feature in operon:
        feature.ignore('')
    result = create_operon_figure(operon, True, None)
    assert result is not None


def test_calculate_adjusted_operon_bounds():
    operon = get_standard_operon()
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


def test_create_operon_figure_with_ignored():
    operon = get_standard_operon()
    fs = FilterSet().must_be_within_n_bp_of_feature('cas2', 10)
    fs.evaluate(operon)
    ax = create_operon_figure(operon, True, None)
    features = _find_plotted_features(ax)
    assert features == set(['cas1', 'cas2', 'cas4 (ignored)'])


def test_create_operon_figure_with_colors():
    operon = get_standard_operon()
    fs = FilterSet().must_be_within_n_bp_of_feature('cas2', 10)
    fs.evaluate(operon)
    gene_colors = {'cas1': 'purple', 'cas2': 'green'}
    ax = create_operon_figure(operon, False, gene_colors)
    features = _find_plotted_features(ax)
    assert features == set(['cas1', 'cas2'])


def test_create_operon_figure():
    operon = get_standard_operon()
    fs = FilterSet().must_be_within_n_bp_of_feature('cas2', 10)
    fs.evaluate(operon)
    ax = create_operon_figure(operon, False, None)
    features = _find_plotted_features(ax)
    assert features == set(['cas1', 'cas2'])
