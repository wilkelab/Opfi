from operon_analyzer import genes, repeat_finder
import pytest


def test_parse_repeats():
    result = '>sequence:14:85:16m\n'
    sequence = "TTCCCCCCTTTTTGGGAATATCGGTATGCTTTTTTTTGATCCCAGAGAGCCCTTTACATAGTTTTTTTTGCATACCGATATTCCCTTCCCCCCTTTTT"
    #                        ++++++++++++++++        oooooooooooooooooooooooo        ----------------
    operon = genes.Operon('Operon', '', 47, 70, [genes.Feature("protein1", (48, 51), '', 1, '', 1e-30, '', 'MQRT', 1230),
                                                 genes.Feature("protein2", (60, 68), '', 1, '', 2e-30, '', 'MGGG', 2401),
                                                 ])
    buffer_size = 37
    start, end, alignment = next(repeat_finder._parse_grf_results(result))
    up_feature, down_feature = repeat_finder._parse_repeats(start, end, alignment, sequence, buffer_size, operon.start, 'perfect')
    assert up_feature.start == 23
    assert up_feature.end == 39
    assert down_feature.start == 79
    assert down_feature.end == 95


@pytest.mark.parametrize('size, expected', [
    (0,   'GGGGGGGGGG'),
    (1,  'AGGGGGGGGGGA'),
    (2, 'TAGGGGGGGGGGAC')
    ])
def test_get_buffered_operon_sequence(size, expected):
    operon = genes.Operon('Operon', '', 11, 20, [genes.Feature("protein", (12, 16), '', 1, '', 1e-30, '', 'MQRT', 1230)])
    sequence = 'AAAAAAAATAGGGGGGGGGGACAAAAAAAA'
    actual = repeat_finder._get_buffered_operon_sequence(operon, sequence, size)
    assert actual == expected


def test_run_grf():
    # The inverted repeats are marked by + and -
    # the operon is marked by o
    sequence = "GGGAATATCGGTATGCTTTTTTTTGATCCCAGAGAGCCCTTTACATAGTTTTTTTTGCATACCGATATTCCC"
    #           ++++++++++++++++        oooooooooooooooooooooooo        ----------------
    perfect, imperfect = repeat_finder._run_grf(sequence, 20, 12)
    assert not imperfect
    assert perfect == '>sequence:1:72:16m'


def test_run_grf_with_surrounding_random_dna():
    # The inverted repeats are marked by + and -
    # the operon is marked by o
    sequence = "TTCCCCCCTTTTTGGGAATATCGGTATGCTTTTTTTTGATCCCAGAGAGCCCTTTACATAGTTTTTTTTGCATACCGATATTCCCTTCCCCCCTTTTT"
    #                        ++++++++++++++++        oooooooooooooooooooooooo        ----------------
    perfect, imperfect = repeat_finder._run_grf(sequence, 20, 12)
    assert not imperfect
    assert perfect == '>sequence:14:85:16m'


def test_parse_repeat_id():
    repeat_id = ">sequence:13:133:54m2D3I"
    start, end, alignment = repeat_finder._parse_repeat_id(repeat_id)
    assert start == 13
    assert end == 133
    assert alignment == "54m2D3I"


@pytest.mark.parametrize('alignment, expected', [
    ('33m', (33, 33)),
    ('17m3D2I1M5m', (25, 26)),
    ('11M11m11M', (33, 33)),
    ('1D1I2M3D7I', (10, 6))])
def test_parse_alignment_size(alignment, expected):
    actual = repeat_finder._parse_alignment_size(alignment)
    assert actual == expected


@pytest.mark.parametrize('upseq, downseq, alignment, expected_up, expected_down', [
     ('AAAGC', 'TTTTT', '5M', 'AAAGC', 'TTTTT'),
     ('AAAGC', 'TTTTT', '5m', 'AAAGC', 'TTTTT'),
     ('AAAGC', 'AAC', '2m2I1m', 'AAAGC', 'AA--C'),
     ('AAAGC', 'AAT', '2m2I1M', 'AAAGC', 'AA--T'),
     ('AAAGC', 'AACTC', '2m2D2I1m', 'AA--AGC', 'AACT--C'),
     ('GGAGC', 'AACTC', '2M2D2I1m', 'GG--AGC', 'AACT--C')])
def test_format_aligned_sequences(upseq, downseq, alignment, expected_up, expected_down):
    actual_up, actual_down = repeat_finder._format_aligned_sequences(upseq, downseq, alignment)
    assert actual_up == expected_up
    assert actual_down == expected_down
