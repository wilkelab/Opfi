from matplotlib.text import Text
from operon_analyzer.genes import Feature, Operon
from operon_analyzer.rules import RuleSet, FilterSet
from operon_analyzer.visualize import calculate_adjusted_operon_bounds, \
                                      create_operon_figure, \
                                      build_image_filename, \
                                      _get_feature_color, \
                                      _make_operon_pairs, \
                                      _plot_clustered_operons
from operon_analyzer.analyze import cluster_operons_by_feature_order
from common import get_standard_operon
import tempfile, os
import pytest


def test_make_operon_pairs():
    features = [Feature('cas1', (0, 20), "", 1, "", 1e-30, "", "MDGYACGYAC", 1234)]

    # perfect overlap
    operon1 = Operon("a", "a.fasta", 0, 100, features)
    other1 = Operon("a", "a.fasta", 0, 100, features)

    # partial overlap
    operon2 = Operon("b", "b.fasta", 200, 300, features)
    other2 = Operon("b", "b.fasta", 250, 350, features)

    # two competing partial overlaps
    operon3 = Operon("c", "c.fasta", 600, 400, features)
    other3good = Operon("c", "c.fasta", 600, 450, features)
    other3bad = Operon("c", "c.fasta", 600, 500, features)

    # no overlap
    operon4 = Operon("d", "d.fasta", 1000, 10000, features)
    other4 = Operon("d", "d.fasta", 20000, 30000, features)

    # other is superset
    operon5 = Operon("e", "e.fasta", 1000, 10000, features)
    other5 = Operon("e", "e.fasta", 400, 12000, features)

    operons = [operon1, operon2, operon3, operon4, operon5]
    others = [other1, other2, other3good, other3bad, other4, other5]

    pairs = _make_operon_pairs(operons, others)
    # Since all the operons must have the same accession and filename for this test to
    # be meaningful, we're distinguishing them by their bounds alone. Hence this janky
    # comparison to expected values:
    expected = {"a": (0, 100),
                "b": (250, 350),
                "c": (600, 450),
                "e": (400, 12000)}
    for operon, candidate in pairs:
        start, end = expected[operon.contig]
        assert start == candidate.start
        assert end == candidate.end
    assert len(pairs) == 4


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


def test_calculate_adjusted_operon_bounds_one_ignored():
    operon = get_standard_operon()
    for feature in operon:
        feature.ignore('')
        break
    result = calculate_adjusted_operon_bounds(operon, False)
    assert result is not None


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
    result = create_operon_figure(operon, False)
    assert result is None


def test_create_operon_figure_all_ignored_plot_ignored():
    operon = get_standard_operon()
    for feature in operon:
        feature.ignore('')
    result = create_operon_figure(operon, True)
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
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)
    fs = FilterSet().must_be_within_n_bp_of_feature('cas2', 10000)
    rs = RuleSet().require('CRISPR array')
    fs.evaluate(operon)
    result = rs.evaluate(operon)
    assert result.is_passing
    ax = create_operon_figure(operon, True)
    features = _find_plotted_features(ax)
    assert features == set(['cas1', 'cas2', 'cas4', 'CRISPR array (2)'])


def test_create_operon_figure_with_ignored():
    operon = get_standard_operon()
    fs = FilterSet().must_be_within_n_bp_of_feature('cas2', 10)
    fs.evaluate(operon)
    ax = create_operon_figure(operon, True)
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
    ax = create_operon_figure(operon, False)
    features = _find_plotted_features(ax)
    assert features == set(['cas1', 'cas2'])


def test_get_feature_color_no_regex():
    feature_colors = {"cas9": "red",
                      "rad51": "green"}
    assert _get_feature_color("cas9", feature_colors) == "red"
    assert _get_feature_color("rad51", feature_colors) == "green"
    assert _get_feature_color("nonexistent", feature_colors) == "blue"


def test_get_feature_color_no_regex_default():
    feature_colors = {"cas9": "red",
                      "rad51": "green",
                      "": "cornflour blue"}
    assert _get_feature_color("cas9", feature_colors) == "red"
    assert _get_feature_color("rad51", feature_colors) == "green"
    assert _get_feature_color("nonexistent", feature_colors) == "cornflour blue"


def test_get_feature_color_regex():
    feature_colors = {"cas9": "red",
                      "rad51": "green"}
    assert _get_feature_color("S. pyogenes Cas9", feature_colors) == "red"
    assert _get_feature_color("Human Rad51", feature_colors) == "green"
    assert _get_feature_color("nonexistent", feature_colors) == "blue"


def test_get_feature_color_regex_default():
    feature_colors = {"cas9": "red",
                      "rad51": "green",
                      "": "cornflour blue"}
    assert _get_feature_color("S. pyogenes Cas9", feature_colors) == "red"
    assert _get_feature_color("Human Rad51", feature_colors) == "green"
    assert _get_feature_color("nonexistent", feature_colors) == "cornflour blue"


def test_get_feature_color_real_regex_default():
    feature_colors = {"cas(9|12a)": "red",
                      "rad51": "green",
                      "": "cornflour blue"}
    assert _get_feature_color("S. pyogenes Cas9, an acceptable, but frankly inferior nuclease", feature_colors) == "red"
    assert _get_feature_color("cas9", feature_colors) == "red"
    assert _get_feature_color("Acidaminococcus sp. BV3L6 Cas12a, CRISPR endonuclease of legend", feature_colors) == "red"
    assert _get_feature_color("Human Rad51", feature_colors) == "green"
    assert _get_feature_color("nonexistent", feature_colors) == "cornflour blue"


@pytest.fixture()
def temporary_directory():
    tmp = tempfile.TemporaryDirectory()
    yield tmp
    tmp.cleanup()


def test_plot_clustered_operons_motif_name_longer_than_sys_limit(temporary_directory):
    # create a dummy operon that has a really long (>255 characters) motif name
    genes = []
    for i in range(100):
        genes.append(Feature('cas1', (i, i+100), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'))
    operons = [Operon('QCDRTU', '/tmp/dna.fasta', 0, 3400, genes)]
    clustered_operons = cluster_operons_by_feature_order(operons)
    _plot_clustered_operons(clustered_operons, image_dir=temporary_directory.name, plot_ignored=True, feature_colors={"cas1": "blue"})
    gene_names = ["cas1"] * 100
    motif_name = "1-" + "-".join(gene_names)
    truncated_motif_name = motif_name[:os.statvfs(temporary_directory.name).f_namemax - 10]
    assert len(os.listdir(temporary_directory.name)) == 1 and os.listdir(temporary_directory.name)[0][:-7] == truncated_motif_name
    