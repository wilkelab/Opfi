from Bio.Seq import Seq
import pytest
from operon_analyzer.piler_parse import RepeatSpacer, BrokenSpacer
from operon_analyzer.spacers import _count_cigar_matches, _align, _build_censored_contig, _fix_broken_spacer


def test_align():
    spacer = "CTGGCCC"
    #         01234567890123456789012
    contig = "AAAAAAAAAAAAAAAACTGCCCAAAAAAAAAAAAAAA"
    _, match_count, contig_target_sequence, contig_start, contig_end, contig_alignment, spacer_alignment, comp_string = _align(spacer, contig)
    assert match_count == 6
    assert contig_alignment == "CT-GCCC"
    assert contig_start == 16
    assert contig_end == 21
    assert spacer_alignment == "CTGGCCC"
    assert comp_string == "|| ||||"


def test_fix_broken_spacer():
    #                 01234567890123456789012345678
    contig_seq = Seq("GGGGGTTTAAAACCCCCCCCCCAAAAAAA")
    # 3 bp repeat ("TTT") followed by 4 bp broken spacer ("AAAA") at position 5
    broken = BrokenSpacer(5, 3, "AAAA")
    fixed = _fix_broken_spacer(broken, 15, contig_seq)
    assert fixed.sequence == "AAAACCCCCCCCCCA"
    assert fixed.position == 5
    assert fixed.repeat_len == broken.repeat_len
    assert fixed.spacer_len == 15


def test_build_censored_contig():
    #             01234567890123456789012345
    contig = Seq("AAAAAAAAAACTGTTGAAAAAAAAAA")
    # 3 bp repeat ("AAA") followed by a 6 bp spacer ("CTGTTG") at position 7
    spacer = RepeatSpacer(7, 3, 6, "CTGTTG")
    actual = _build_censored_contig(spacer, contig)
    expected = "AAAAAAAAAANNNNNNAAAAAAAAAA"
    assert actual == expected
    assert len(actual) == len(contig)


@pytest.mark.parametrize('string,expected', [
    ('10=', 10),
    ('3=', 3),
    ('124=', 124),
    ('17=1I4=', 21),
    ('3X14=', 14),
    ('14=3X', 14),
    ('4X1D4=10I7=', 11)])
def test_count_cigar_matches(string, expected):
    actual = _count_cigar_matches(string)
    assert actual == expected
