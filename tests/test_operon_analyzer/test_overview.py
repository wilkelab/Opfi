from operon_analyzer.overview import _count_results, _extract_results, load_counts


def test_load_counts():
    csv_text = [
            'ABC,1,299,fail,exclude:cas3,max-distance-to-anything:transposase-500,min-distance-to-anything:transposase-1',
            'DEF,2,383,fail,require:transposase,max-distance-to-anything:transposase-500,min-distance-to-anything:transposase-1',
            'GHI,345,678,pass',
            'JKL,230,1000,fail,exclude:cas3',
            'MNO,1200,1400,fail,require:transposase',
            'PQR,200,1300,fail,exclude:cas3,max-distance-to-anything:transposase-500,min-distance-to-anything:transposase-1',
            'STU,100,12012,fail,max-distance-to-anything:transposase-500',
            'TUV,300,900,fail,exclude:cas3']
    unique_rule_violated, failed_rule_occurrences, rule_failure_counts = load_counts(csv_text)
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


def test_extract_results():
    lines = [["YFEFWEF323", 390, 1231241, ["pass"]],
             ["FYYFF233", 44858, 12301231, ["fail", "require:cas2", "exclude:cas3"]]]
    actual = list(_extract_results(lines))
    expected = [["pass"], ["fail", "require:cas2", "exclude:cas3"]]
    assert actual == expected


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
