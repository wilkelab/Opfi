from operon_analyzer.genes import Feature, Operon
from operon_analyzer import parse
from Bio.Seq import Seq


def test_serialize_operon_roundtrip():
    line = "forward,473846..494916,cas5,485260..488403,lcl|485260|488403|1|1,1,UniRef50_D3I2J9,7.28e-113,UniRef50_D3I2J9 cas5 Putative CRISPR-associated protein Csc1 n=138 RepID=D3I2J9_9BACT,MNLTLKTLLALNLTLI,374,959,1037,29.219,303,617,498,30,117,48.02,93,forward.fa\n"
    operon = list(parse.load_operons([line]))[0]
    actual = operon.as_str()
    assert actual == line


def test_feature_length():
    f = Feature('gene', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    assert len(f) == 100


def test_operon_feature_region_sequence():
    f1 = Feature('gene1', (5, 10), '', 1, '', 4e-12, '', 'MG', 1234.9)
    f2 = Feature('gene2', (12, 20), '', 1, '', 4e-12, '', 'MTE', 1234.9)
    f3 = Feature('gene3', (25, 30), '', 1, '', 4e-12, '', 'MR', 1234.9)
    op = Operon('asdf', '/tmp/dna.fasta', 0, 999, [f1, f2, f3])
    # 01234567890123456789012345678901
    # AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT
    #      +++++  ++++++++    ++++++
    op.set_sequence(Seq("AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT"))
    expected = "AAACCCCCCCCGGGGGGGGTTTTTT"
    assert op.feature_region_sequence == expected


def test_get_unique():
    f1 = Feature('gene1', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f2 = Feature('gene2', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f3 = Feature('gene3', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f4 = Feature('gene3', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    op = Operon('asdf', '/tmp/dna.fasta', 0, 1000, [f1, f2, f3, f4])
    funique = op.get_unique('gene1')
    assert funique is not None


def test_get_unique_not_unique():
    f1 = Feature('gene1', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f2 = Feature('gene2', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f3 = Feature('gene3', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f4 = Feature('gene3', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    op = Operon('asdf', '/tmp/dna.fasta', 0, 1000, [f1, f2, f3, f4])
    fnotunique = op.get_unique('gene3')
    assert fnotunique is None


def test_feature_repr():
    f = Feature('gene1', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    label = str(f)
    assert label == '<Feature gene1 1..100>'


def test_redundant_features():
    f1 = Feature('gene1', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f2 = Feature('gene1', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f3 = Feature('gene2', (1, 100), '', 1, '', 4e-12, '', 'MQQHT', 1234.9)
    f4 = Feature('gene3', (1, 100), '', 1, '', 4e-12, '', 'MRTKKP', 1234.9)
    s = set([f1, f2, f3, f4])
    assert len(s) == 3


def test_identical_crispr_arrays_redundant():
    f1 = Feature('CRISPR array', (1300, 1400), '', None, '', None, 'Copies: 2, Repeat: 61, Spacer: 59', 'CTCAAACTCTCACTCTGGCTCAATGAGTTAGACAAGCTCTCACTCTGACTCAAGGAATTAC')
    f2 = Feature('CRISPR array', (1300, 1400), '', None, '', None, 'Copies: 2, Repeat: 61, Spacer: 59', 'CTCAAACTCTCACTCTGGCTCAATGAGTTAGACAAGCTCTCACTCTGACTCAAGGAATTAC')
    s = set([f1, f2])
    assert len(s) == 1


def test_nonidentical_crispr_arrays_nonredundant():
    f1 = Feature('CRISPR array', (1305, 1405), '', None, '', None, 'Copies: 2, Repeat: 61, Spacer: 59', 'CTCAAACTCTCACTCTGGCTCAATGAGTTAGACAAGCTCTCACTCTGACTCAAGGAATTAC')
    f2 = Feature('CRISPR array', (1300, 1400), '', None, '', None, 'Copies: 2, Repeat: 61, Spacer: 59', 'CTCAAACTCTCACTCTGGCTCAATGAGTTAGACAAGCTCTCACTCTGACTCAAGGAATTAC')
    s = set([f1, f2])
    assert len(s) == 2


def test_nonidentical_crispr_arrays_nonredundant2():
    f1 = Feature('CRISPR array', (1300, 1400), '', None, '', None, 'Copies: 2, Repeat: 61, Spacer: 59', 'CTCAAACTCTCACTCTGGCTCAATGAGTTAGACAAGCTCTCACTCTGACTCAAGGAATTAC')
    f2 = Feature('CRISPR array', (1300, 1400), '', None, '', None, 'Copies: 2, Repeat: 61, Spacer: 59', 'GGGGGGCTCTCACTCTGGCTCAATGAGTTAGACAAGCTCTCACTCTGACTCAAGGAATTAC')
    s = set([f1, f2])
    assert len(s) == 2


def test_redundant_operons():
    # Operons are identical except for different source filename
    f1 = Feature('gene1', (1, 100), '', 1, 'HFEFIFJF3.11', 4e-12, 'Methyltransferase', 'MGWRN', 1234.9)
    f2 = Feature('gene1', (1, 100), '', 1, 'HFEFIFJF3.11', 4e-12, 'Methyltransferase', 'MGWRN', 1234.9)
    f3 = Feature('gene2', (200, 300), '', 1, 'AMGFEFJEFI1.2', 7e-12, 'DNA polymerase', 'MQQHTR', 1098.9)
    f4 = Feature('gene2', (200, 300), '', 1, 'AMGFEFJEFI1.2', 7e-12, 'DNA polymerase', 'MQQHTR', 1098.9)
    f5 = Feature('CRISPR array', (1300, 1400), '', None, '', None, 'Copies: 2, Repeat: 61, Spacer: 59', 'CTCAAACTCTCACTCTGGCTCAATGAGTTAGACAAGCTCTCACTCTGACTCAAGGAATTAC')
    f6 = Feature('CRISPR array', (1300, 1400), '', None, '', None, 'Copies: 2, Repeat: 61, Spacer: 59', 'CTCAAACTCTCACTCTGGCTCAATGAGTTAGACAAGCTCTCACTCTGACTCAAGGAATTAC')
    op1 = Operon('AA3.1', '/path/to/operon.fa.gz', 450, 10000, [f1, f3, f5])
    op2 = Operon('AA3.1', '/alternate/path/to/operon.fa.gz', 450, 10000, [f2, f4, f6])
    s = set([op1, op2])
    assert len(s) == 1


def test_nonredundant_operons():
    # Operons are identical except one has one less Feature than the other
    f1 = Feature('gene1', (1, 100), '', 1, 'HFEFIFJF3.11', 4e-12, 'Methyltransferase', 'MGWRN', 1234.9)
    f3 = Feature('gene1', (1, 100), '', 1, 'HFEFIFJF3.11', 4e-12, 'Methyltransferase', 'MGWRN', 1234.9)
    f2 = Feature('gene2', (200, 300), '', 1, 'AMGFEFJEFI1.2', 7e-12, 'DNA polymerase', 'MQQHTR', 1098.9)
    op1 = Operon('AA3.1', '/path/to/operon.fa.gz', 450, 10000, [f1, f2])
    op2 = Operon('AA3.1', '/path/to/operon.fa.gz', 450, 10000, [f3])
    s = set([op1, op2])
    assert len(s) == 2


def test_nonredundant_operons2():
    # Operons vary only by accession ID
    f1 = Feature('gene1', (1, 100), '', 1, 'HFEFIFJF3.11', 4e-12, 'Methyltransferase', 'MGWRN', 1234.9)
    f3 = Feature('gene1', (1, 100), '', 1, 'HFEFIFJF3.11', 4e-12, 'Methyltransferase', 'MGWRN', 1234.9)
    f2 = Feature('gene2', (200, 300), '', 1, 'AMGFEFJEFI1.2', 7e-12, 'DNA polymerase', 'MQQHTR', 1098.9)
    f4 = Feature('gene2', (200, 300), '', 1, 'AMGFEFJEFI1.2', 7e-12, 'DNA polymerase', 'MQQHTR', 1098.9)
    op1 = Operon('AA3.1', '/path/to/operon.fa.gz', 450, 10000, [f1, f2])
    op2 = Operon('ZZ3.1', '/path/to/operon.fa.gz', 450, 10000, [f3, f4])
    s = set([op1, op2])
    assert len(s) == 2
