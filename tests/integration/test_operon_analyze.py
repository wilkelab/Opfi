import tempfile
import os
import shutil
from operon_analyzer.analyze import analyze
from operon_analyzer.rules import RuleSet
from operon_analyzer.visualize import build_operon_dictionary, load_analyzed_operons, plot_operons


def test_analyze(capsys):
    """ Just serves to check that `analyze()` produces output """
    rs = RuleSet().require('transposase') \
                  .exclude('cas3') \
                  .max_distance_to_anything('transposase', 500) \
                  .min_distance_to_anything('transposase', 1)

    with open('tests/integration/integration_data/operon_analyzer/transposases.csv') as f:
        analyze(f, rs)
        captured = capsys.readouterr()
        stdout = captured.out
        assert stdout.startswith("#")
        assert stdout.count("pass") == 1


def test_visualize_passes():
    pass_count = visualize('pass')
    assert pass_count == 2


def test_visualize_failures():
    fail_count = visualize('fail')
    assert fail_count == 4


def test_visualize_all():
    count = visualize('')
    assert count == 6


def test_visualize_none():
    count = visualize('nonexistent-condition')
    assert count == 0


def visualize(condition: str):
    """
    Creates PNGs of the operons matching the given condition. The idea here is that
    there are four failing and two passing operons. We make PNGs of operons whose
    analysis result field starts with the given condition (which is either "pass" or "fail").
    Tests that use this function just determine whether the expected number of PNGs were made,
    not whether they are correct.
    """
    analysis_csv = 'tests/integration/integration_data/operon_analyzer/analysis.csv'
    pipeline_csv = 'tests/integration/integration_data/operon_analyzer/pipeline.csv'

    # We make a temporary directory to store the PNGs
    tempdir = tempfile.mkdtemp()
    try:
        good_operons = []
        with open(pipeline_csv) as f:
            operons = build_operon_dictionary(f)
        with open(analysis_csv) as f:
            for contig, start, end, result in load_analyzed_operons(f):
                if not result.startswith(condition):
                    continue
                op = operons.get((contig, start, end))
                if op is None:
                    continue
                good_operons.append(op)
        plot_operons(good_operons, tempdir)
        files = os.listdir(tempdir)
        count = len([f for f in files if f.endswith(".png")])
    except Exception as e:
        raise e
    finally:
        # clean up the directory and any PNGs that were made
        shutil.rmtree(tempdir)
    return count
