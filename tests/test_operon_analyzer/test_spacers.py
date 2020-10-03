from operon_analyzer.spacers import _count_cigar_matches
import pytest


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
