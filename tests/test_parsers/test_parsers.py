import pytest
from crisposon.parsers import (
    parse_diamond,
    parse_mmseqs
)

def ret_mmseqs_fields():

    fields = ["query","target","evalue","qseq",
                "qheader","theader","qcov","tset"]

    return fields

def ret_diamond_fields():

    fields = ["qseqid", "sseqid", "full_qseq", "evalue", "stitle"]
    return fields

def test_mmseqs_empty_tsv():
    """Test that parse_mmseqs returns an empty dict 
    if the tsv file is empty.
    """

    fields = ret_mmseqs_fields()
    result = parse_mmseqs("tests/test_parsers/tsv/empty.tsv", 
                            "empty", fields)

    assert isinstance(result, dict)
    assert len(result) == 0

def test_mmseqs():
    """Test parse_mmseqs on a tsv file with multiple hits
    for multiple queries.
    """

    fields = ret_mmseqs_fields()
    result = parse_mmseqs("tests/test_parsers/tsv/mmseqs.tsv",
                            "mmseqs", fields)

    assert len(result) == 2
    
    evalues = [float(hit["Hit_e-val"]) for hit in result.values()]
    assert 8.654E-87 and 1.992E-27 in evalues

def test_diamond_empty_tsv():
    """Test that parse_diamond returns an empty dict 
    if the output tsv file is empty.
    """

    fields = ret_diamond_fields()
    result = parse_diamond("tests/test_parsers/tsv/empty.tsv",
                            "empty", fields)

    assert isinstance(result, dict)
    assert len(result) == 0

def test_diamond_no_des():
    """Test parse_diamond on an output tsv file generated
    with a reference database with that doesn't have
    headers.
    """

    fields = ret_diamond_fields()
    result = parse_diamond("tests/test_parsers/tsv/diamond_no_des.tsv",
                            "diamond", fields)
    
    for hit in result.values():
        assert hit["Hit_accession"] == ''
        assert hit["Hit_description"] == ''

def test_diamond():
    """Test parse_diamond on a tsv file with multiple hits
    for multiple queries.
    """

    fields = ret_diamond_fields()
    result = parse_diamond("tests/test_parsers/tsv/diamond_des.tsv",
                             "diamond", fields)

    assert len(result) == 2
    
    evalues = [float(hit["Hit_e-val"]) for hit in result.values()]
    assert 5.8e-155 and 1.3e-119 in evalues
