from operon_analyzer.analyze import analyze
from operon_analyzer.rules import RuleSet


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
