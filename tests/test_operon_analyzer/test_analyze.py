from operon_analyzer.genes import Feature, Operon
from operon_analyzer import analyze


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
    actual = analyze._get_sorted_feature_names(operon)
    expected = ('cas2', 'cas4', 'cas1')
    assert actual == expected


def test_sort_feature_names3():
    genes = [
            Feature('cas2', (4000, 1200), 'lcl|12|400|1|-1', 1, 'ACACEHFEF', 4e-19, 'a good gene', 'MCGYVER'),
            Feature('cas1', (410, 600), 'lcl|410|600|1|-1', 1, 'FGEYFWCE', 2e-5, 'a good gene', 'MGFRERAR'),
            Feature('cas4', (620, 2200), 'lcl|620|1200|1|-1', 1, 'NFBEWFUWEF', 6e-13, 'a good gene', 'MLAWPVTLE'),
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

