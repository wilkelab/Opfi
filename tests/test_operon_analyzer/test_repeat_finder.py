from Bio.Seq import Seq
from operon_analyzer import genes, repeat_finder
import pytest


def test_find_inverted_repeats():
    contig_sequence = Seq("TTCCCCCCTTTTTGGGAATATCGGTATGCTTTTTTTTGATCCCAGAGAGCCCTTTACATAGTTTTTTTTGCATACCGATATTCCCTTCCCCCCTTTTT")
    operon = genes.Operon('Operon', '', 37, 61, [genes.Feature("protein1", (48, 51), '', 1, '', 1e-30, '', 'M', 1230),
                                                 genes.Feature("protein2", (54, 61), '', 1, '', 2e-30, '', 'MG', 2401)])
    repeat_finder._find_inverted_repeats(operon, contig_sequence, 40, 12)
    irs = [feat for feat in operon if feat.name.startswith("IR")]
    assert len(irs) == 2
    assert irs[0].sequence == 'GGGAATATCGGTATGC'
    assert irs[1].sequence == Seq('GCATACCGATATTCCC').reverse_complement()


def test_buffered_sequence():
    #                                10        20        30        40        50        60        70        80        90
    #                      01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567
    #                              0         10        20        30        40        50        60        70        80        90
    contig_sequence = Seq("TTCCCCCCTTTTTGGGAATATCGGTATGCTTTTTTTTGATCCCAGAGAGCCCTTTACATAGTTTTTTTTGCATACCGATATTCCCTTCCCCCCTTTTT")
    #                                   ++++++++++++++++        oooooooooooooooooooooooo        ----------------
    operon = genes.Operon('Operon', '', 37, 61, [genes.Feature("protein1", (48, 51), '', 1, '', 1e-30, '', 'M', 1230),
                                                 genes.Feature("protein2", (54, 61), '', 1, '', 2e-30, '', 'MG', 2401),
                                                 ])
    bufseq = repeat_finder.BufferedSequence(operon, contig_sequence, 40)
    expected_seq = "TTTTTGGGAATATCGGTATGCTTTTTTTTGATCCCAGAGAGCCCTTTACATAGTTTTTTTTGCATACCGATATTCCCTTCCCCCCTTTTT"
    assert bufseq.start == 8
    assert bufseq.end == 98
    assert bufseq.sequence == expected_seq
    assert bufseq[13:29] == "CGGTATGCTTTTTTTT"
    assert bufseq[69:85] == "ATATTCCCTTCCCCCC"
    assert bufseq.operon_length == 14


def test_parse_repeats():
    result = '>sequence:14:85:16m\n'
    #                                10        20        30        40        50        60        70        80        90
    #                      01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567
    contig_sequence = Seq("TTCCCCCCTTTTTGGGAATATCGGTATGCTTTTTTTTGATCCCAGAGAGCCCTTTACATAGTTTTTTTTGCATACCGATATTCCCTTCCCCCCTTTTT")
    #                               ++++++++++++++++        oooooooooooooooooooooooo        ----------------
    operon = genes.Operon('Operon', '', 37, 61, [genes.Feature("protein1", (48, 51), '', 1, '', 1e-30, '', 'M', 1230),
                                                 genes.Feature("protein2", (54, 61), '', 1, '', 2e-30, '', 'MG', 2401),
                                                 ])
    buffer_size = 60
    bufseq = repeat_finder.BufferedSequence(operon, contig_sequence, buffer_size)
    result = next(repeat_finder._parse_grf_results(result))
    up_feature, down_feature = repeat_finder._parse_repeats(result, bufseq, 77)
    assert up_feature.start == 13
    assert up_feature.end == 29
    assert '77' in up_feature.name
    assert down_feature.start == 69
    assert down_feature.end == 85
    assert '77' in down_feature.name


def test_run_grf():
    # The inverted repeats are marked by + and -
    # the operon is marked by o
    operon = genes.Operon('Operon', '', 24, 48, [genes.Feature("protein1", (24, 27), '', 1, '', 1e-30, '', 'M', 1230),
                                                 genes.Feature("protein2", (39, 48), '', 1, '', 2e-30, '', 'MG', 2401),
                                                 ])
    #                                                    0         10        20        30        40        50        60        70
    bufseq = repeat_finder.BufferedSequence(operon, Seq("GGGAATATCGGTATGCTTTTTTTTGATCCCAGAGAGCCCTTTACATAGTTTTTTTTGCATACCGATATTCCC"), 30)
    #                                                    ++++++++++++++++        oooooooooooooooooooooooo        ----------------
    perfect, imperfect = repeat_finder._run_grf(bufseq, 12)
    assert not imperfect
    assert perfect == '>sequence:1:72:16m'


def test_run_grf_with_surrounding_random_dna():
    # The inverted repeats are marked by + and -
    # the operon is marked by o
    #                      0         10        20        30        40        50        60        70        80        90     
    contig_sequence = Seq("TTCCCCCCTTTTTGGGAATATCGGTATGCTTTTTTTTGATCCCAGAGAGCCCTTTACATAGTTTTTTTTGCATACCGATATTCCCTTCCCCCCTTTTT")
    #                                   ++++++++++++++++        oooooooooooooooooooooooo        ----------------
    operon = genes.Operon('op123', '/path/to/fa.gz', 37, 61, [genes.Feature('gene', (37, 61), '', 1, '', 1e-30, '', 'MC', 123)])
    buffered_sequence = repeat_finder.BufferedSequence(operon, contig_sequence, 37)
    perfect, imperfect = repeat_finder._run_grf(buffered_sequence, 12)
    assert not imperfect
    assert perfect == '>sequence:14:85:16m'


def test_run_grf_with_surrounding_random_dna_shortened_buffer():
    # The inverted repeats are marked by + and -
    # the operon is marked by o
    #                      0         10        20        30        40        50        60        70        80        90     
    contig_sequence = Seq("TTCCCCCCTTTTTGGGAATATCGGTATGCTTTTTTTTGATCCCAGAGAGCCCTTTACATAGTTTTTTTTGCATACCGATATTCCCTTCCCCCCTTTTT")
    #                                   ++++++++++++++++        oooooooooooooooooooooooo        ----------------
    operon = genes.Operon('op123', '/path/to/fa.gz', 37, 61, [genes.Feature('gene', (37, 61), '', 1, '', 1e-30, '', 'MC', 123)])
    buffered_sequence = repeat_finder.BufferedSequence(operon, contig_sequence, 35)
    perfect, imperfect = repeat_finder._run_grf(buffered_sequence, 12)
    assert not imperfect
    assert perfect == '>sequence:12:83:16m'


def test_parse_repeat_id():
    repeat_id = ">sequence:13:133:54m2D3I"
    start, end, alignment = repeat_finder._parse_repeat_id(repeat_id)
    assert start == 12
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
