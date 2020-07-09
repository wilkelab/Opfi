import pytest
from gene_finder.parsers import (
    parse_search_output
)

def test_empty_tsv():
    """
    Test that parser correctly handles an empty tsv.
    """

    result = parse_search_output("tests/test_gene_finder/tsv/empty.tsv", 
                            "empty", "blast")

    assert isinstance(result, dict)
    assert len(result) == 0


def test_mmseqs():
    """Test parse_mmseqs on a tsv file with multiple hits
    for multiple queries.
    """

    result = parse_search_output("tests/test_gene_finder/tsv/mmseqs_parse_test.tsv",
                            "mmseqs", "mmseqs")

    assert len(result) == 3
    
    evalues = [float(hit["Hit_e-val"]) for hit in result.values()]
    assert 5.991E-04 and 3.134E-05 and 1.992E-27 in evalues


def test_diamond():
    """Test parse_diamond on a tsv file with multiple hits
    for multiple queries.
    """

    result = parse_search_output("tests/test_gene_finder/tsv/diamond_parse_test.tsv",
                             "diamond", "diamond")

    assert len(result) == 2
    
    evalues = [float(hit["Hit_e-val"]) for hit in result.values()]
    assert 1.2e-22 and 1.3e-119 in evalues
