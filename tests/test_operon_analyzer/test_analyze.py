from operon_analyzer.genes import Feature, Operon
from operon_analyzer import analyze
from Bio.Seq import Seq


def test_operons_approximately_equal():
    genes1 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon1 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes1)
    genes2 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon2 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes2)
    assert analyze._operons_are_approximately_equal(operon1, operon2)


def test_operons_approximately_equal_reversed():
    genes1 = [Feature('cas2', (15, 25), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
              Feature('cas1', (55, 65), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             ]
    operon1 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes1)
    genes2 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon2 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes2)
    assert analyze._operons_are_approximately_equal(operon1, operon2)


def test_operons_approximately_equal_different_feature_position():
    genes1 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon1 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes1)
    genes2 = [Feature('cas1', (20, 30), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon2 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes2)
    assert not analyze._operons_are_approximately_equal(operon1, operon2)


def test_operons_approximately_equal_different_feature_seq1():
    genes1 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERARQQQ')]
    operon1 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes1)
    genes2 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon2 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes2)
    assert not analyze._operons_are_approximately_equal(operon1, operon2)


def test_operons_approximately_equal_different_feature_seq2():
    genes1 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVERQQQ'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon1 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes1)
    genes2 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon2 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes2)
    assert not analyze._operons_are_approximately_equal(operon1, operon2)


def test_group_redundant_operons():
    operons = []
    # The first three operons should all group together
    genes = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon = Operon('A', '/tmp/dna.fasta', 0, 13400, genes)
    operon.set_sequence(Seq("A" * 16 + "C" * 16 + "C" * 16 + "G" * 16))
    operons.append(operon)

    genes = [Feature('cas1', (12, 22), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (52, 62), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon = Operon('B', '/tmp/dna.fasta', 0, 13400, genes)
    operon.set_sequence(Seq("GG" + "A" * 16 + "C" * 16 + "C" * 16 + "G" * 16))
    operons.append(operon)

    genes = [Feature('cas1', (55, 45), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (15, 5), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon = Operon('C', '/tmp/dna.fasta', 0, 13400, genes)
    operon.set_sequence(Seq("C" * 16 + "G" * 16 + "G" * 16 + "T" * 16))
    operons.append(operon)

    # This operon has Features in the exact same location as Operon A, but a different nucleotide sequence
    genes = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon = Operon('D', '/tmp/dna.fasta', 0, 13400, genes)
    operon.set_sequence(Seq("T" * 16 + "T" * 16 + "T" * 16 + "T" * 16))
    operons.append(operon)

    # This operon has the same nucleotide sequence as Operon A, but Features in a different location
    genes = [Feature('cas1', (5, 15), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (45, 55), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon = Operon('E', '/tmp/dna.fasta', 0, 13400, genes)
    operon.set_sequence(Seq("A" * 16 + "C" * 16 + "C" * 16 + "G" * 16))
    operons.append(operon)

    unique = analyze.group_similar_operons(operons, load_sequences=False)
    actual_ids = [operon.contig for operon in unique]
    expected_ids = ["A", "D", "E"]
    assert actual_ids == expected_ids


def test_sort_feature_names():
    genes = [
            Feature('cas1', (1210, 4000), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 13400, genes)
    actual = analyze._get_sorted_feature_names(operon)
    expected = ('cas2', 'cas4', 'cas1')
    assert actual == expected


def test_sort_feature_names2():
    genes = [
            Feature('cas1', (4000, 1200), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 2200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 13400, genes)
def test_operons_approximately_equal_different_feature_position():
    genes1 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon1 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes1)
    genes2 = [Feature('cas1', (20, 30), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon2 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes2)
    assert not analyze._operons_are_approximately_equal(operon1, operon2)


def test_operons_approximately_equal_different_feature_seq1():
    genes1 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERARQQQ')]
    operon1 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes1)
    genes2 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon2 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes2)
    assert not analyze._operons_are_approximately_equal(operon1, operon2)


def test_operons_approximately_equal_different_feature_seq2():
    genes1 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVERQQQ'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon1 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes1)
    genes2 = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon2 = Operon('A', '/tmp/dna.fasta', 0, 13400, genes2)
    assert not analyze._operons_are_approximately_equal(operon1, operon2)


def test_group_redundant_operons():
    operons = []
    # The first three operons should all group together
    genes = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon = Operon('A', '/tmp/dna.fasta', 0, 13400, genes)
    operon.set_sequence(Seq("A" * 16 + "C" * 16 + "C" * 16 + "G" * 16))
    operons.append(operon)

    genes = [Feature('cas1', (12, 22), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (52, 62), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon = Operon('B', '/tmp/dna.fasta', 0, 13400, genes)
    operon.set_sequence(Seq("GG" + "A" * 16 + "C" * 16 + "C" * 16 + "G" * 16))
    operons.append(operon)

    genes = [Feature('cas1', (45, 55), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (5, 15), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon = Operon('C', '/tmp/dna.fasta', 0, 13400, genes)
    operon.set_sequence(Seq("C" * 16 + "G" * 16 + "G" * 16 + "T" * 16))
    operons.append(operon)

    # This operon has Features in the exact same location as Operon A, but a different nucleotide sequence
    genes = [Feature('cas1', (10, 20), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (50, 60), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon = Operon('D', '/tmp/dna.fasta', 0, 13400, genes)
    operon.set_sequence(Seq("T" * 16 + "T" * 16 + "T" * 16 + "T" * 16))
    operons.append(operon)

    # This operon has the same nucleotide sequence as Operon A, but Features in a different location
    genes = [Feature('cas1', (5, 15), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
             Feature('cas2', (45, 55), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR')]
    operon = Operon('E', '/tmp/dna.fasta', 0, 13400, genes)
    operon.set_sequence(Seq("A" * 16 + "C" * 16 + "C" * 16 + "G" * 16))
    operons.append(operon)

    unique = analyze.group_similar_operons(operons, load_sequences=False)
    actual_ids = [operon.contig for operon in unique]
    expected_ids = ["A", "D", "E"]
    assert actual_ids == expected_ids


def test_sort_feature_names():
    genes = [
            Feature('cas1', (1210, 4000), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 1200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 13400, genes)
    actual = analyze._get_sorted_feature_names(operon)
    expected = ('cas2', 'cas4', 'cas1')
    assert actual == expected


def test_sort_feature_names2():
    genes = [
            Feature('cas1', (1200, 4000), 'lcl|4000|1200|1|-1', -1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (410, 600), 'lcl|410|600|1|1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 2200), 'lcl|620|2200|1|1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 13400, genes)
    actual = analyze._get_sorted_feature_names(operon)
    expected = ('cas2', 'cas4', 'cas1')
    assert actual == expected


def test_sort_feature_names3():
    genes = [
            Feature('cas2', (1200, 4000), 'lcl|4000|1200|1|-1', -1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas1', (410, 600), 'lcl|410|600|1|1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 2200), 'lcl|620|2200|1|1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon = Operon('QCDRTU', '/tmp/dna.fasta', 0, 13400, genes)
    actual = analyze._get_sorted_feature_names(operon)
    expected = ('cas1', 'cas4', 'cas2')
    assert actual == expected


def test_cluster_operons_by_feature_order():
    genes = [
            Feature('cas1', (1, 100), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (200, 300), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (400, 500), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon1 = Operon('QCDRTU', '/tmp/dna.fasta', 0, 13400, genes)
    genes = [
            Feature('cas4', (1, 100), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (200, 300), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas1', (400, 500), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon2 = Operon('QCDRTU', '/tmp/dna.fasta', 0, 13400, genes)
    genes = [
            Feature('cas1', (1, 100), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (200, 300), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (400, 500), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon3 = Operon('QCDRTU', '/tmp/dna.fasta', 0, 13400, genes)
    genes = [
            Feature('cas1', (1, 100), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas4', (400, 500), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon4 = Operon('QCDRTU', '/tmp/dna.fasta', 0, 13400, genes)
    genes = [
            Feature('transposase', (1, 100), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas2', (200, 300), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (400, 500), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
            ]
    operon5 = Operon('QCDRTU', '/tmp/dna.fasta', 0, 13400, genes)
    operons = [operon1, operon2, operon3, operon4, operon5]
    classified = {name: len(ops) for name, ops in analyze.cluster_operons_by_feature_order(operons).items()}
    expected = {('cas1', 'cas2', 'cas4'): 3, ('cas1', 'cas4'): 1, ('transposase', 'cas2', 'cas4'): 1}
    assert classified == expected

