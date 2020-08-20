from operon_analyzer import reannotation, genes
from common import get_standard_operon


def test_reannotate_operon():
    operon = get_standard_operon()

    f1 = genes.Feature('CRISPR Cas1', (12, 400), '', 1, '', 4e-20, 'awesome protein', 'MCGYVERS')
    f2 = genes.Feature('DNA polymerase', (410, 600), '', 1, '', 4e-20, 'awesome protein', 'MRRQTHPR')
    f3 = genes.Feature('CRISPR Cas4', (620, 1200), '', 1, '', 4e-20, 'awesome protein', 'MLCDWALAP')
    reannotated_operon = genes.Operon('operon', '/path/to/contig.fa.gz', 0, 3400, [f1, f2, f3])
    actual = reannotation.count_reannotations(operon, reannotated_operon)
    expected = {('cas1', 1): {'CRISPR Cas1': 1},
                ('cas2', 1): {'DNA polymerase': 1},
                ('cas4', 1): {'CRISPR Cas4': 1}}
    assert actual == expected


def test_reannotate_operon_multi_reannotations():
    operon = get_standard_operon()

    f1 = genes.Feature('CRISPR Cas1', (12, 400), '', 1, '', 4e-20, 'awesome protein', 'MCGYVERS')
    f2 = genes.Feature('DNA polymerase', (410, 500), '', 1, '', 4e-20, 'awesome protein', 'MRRQTHPR')
    f2b = genes.Feature('Methyltransferase', (500, 600), '', 1, '', 4e-20, 'awesome protein', 'MRRQTHPR')
    f3 = genes.Feature('CRISPR Cas4', (620, 1200), '', 1, '', 4e-20, 'awesome protein', 'MLCDWALAP')
    reannotated_operon = genes.Operon('operon', '/path/to/contig.fa.gz', 0, 3400, [f1, f2, f2b, f3])
    actual = reannotation.count_reannotations(operon, reannotated_operon)
    expected = {('cas1', 1): {'CRISPR Cas1': 1},
                ('cas2', 1): {'DNA polymerase': 1, 'Methyltransferase': 1},
                ('cas4', 1): {'CRISPR Cas4': 1}}
    assert actual == expected


def test_reannotate_operon_multi_annotations():
    operon = get_standard_operon()

    f1 = genes.Feature('CRISPR Cas1', (12, 400), '', 1, '', 4e-20, 'awesome protein', 'MCGYVERS')
    f2 = genes.Feature('DNA polymerase', (410, 1200), '', 1, '', 4e-20, 'awesome protein', 'MRRQTHPR')
    reannotated_operon = genes.Operon('operon', '/path/to/contig.fa.gz', 0, 3400, [f1, f2])
    actual = reannotation.count_reannotations(operon, reannotated_operon)
    expected = {('cas1', 1): {'CRISPR Cas1': 1},
                ('cas2', 1): {'DNA polymerase': 1},
                ('cas4', 1): {'DNA polymerase': 1}}
    assert actual == expected


def test_count_cluster_reannotations():
    operon = get_standard_operon()
    operon.contig = 'operon1'
    operon2 = get_standard_operon()
    operon2.contig = 'operon2'

    f1 = genes.Feature('CRISPR Cas1', (12, 400), '', 1, '', 4e-20, 'awesome protein', 'MCGYVERS')
    f2 = genes.Feature('DNA polymerase', (410, 600), '', 1, '', 4e-20, 'awesome protein', 'MRRQTHPR')
    f3 = genes.Feature('CRISPR Cas4', (620, 1200), '', 1, '', 4e-20, 'awesome protein', 'MLCDWALAP')
    reannotated_operon = genes.Operon('operon1', '/path/to/contig.fa.gz', 0, 3400, [f1, f2, f3])

    f1 = genes.Feature('CRISPR Cas1', (12, 400), '', 1, '', 4e-20, 'awesome protein', 'MCGYVERS')
    f2 = genes.Feature('RNA polymerase', (410, 601), '', 1, '', 4e-20, 'awesome protein', 'MRRQTHPR')
    f3 = genes.Feature('CRISPR Cas4 endonuclease', (620, 1205), '', 1, '', 4e-20, 'awesome protein', 'MLCDWALAP')
    reannotated_operon2 = genes.Operon('operon2', '/path/to/contig.fa.gz', 0, 3400, [f1, f2, f3])

    actual = reannotation.count_cluster_reannotations([operon, operon2], {'operon1': reannotated_operon, 'operon2': reannotated_operon2})
    expected = {('cas1', 1): {'CRISPR Cas1': 2},
                ('cas2', 1): {'DNA polymerase': 1, 'RNA polymerase': 1},
                ('cas4', 1): {'CRISPR Cas4': 1, 'CRISPR Cas4 endonuclease': 1}}
    assert actual == expected


def test_count_cluster_reannotations_double_gene():
    operon = get_standard_operon()
    extra_cas1 = genes.Feature('cas1', (1210, 1300), '', 1, '', 4e-20, 'awesome protein', 'MCGYVERS')
    operon._features.append(extra_cas1)
    operon.contig = 'operon1'
    operon2 = get_standard_operon()
    operon2.contig = 'operon2'
    operon2._features.append(extra_cas1)

    f1 = genes.Feature('CRISPR Cas1', (12, 400), '', 1, '', 4e-20, 'awesome protein', 'MCGYVERS')
    f2 = genes.Feature('DNA polymerase', (410, 600), '', 1, '', 4e-20, 'awesome protein', 'MRRQTHPR')
    f3 = genes.Feature('CRISPR Cas4', (620, 1200), '', 1, '', 4e-20, 'awesome protein', 'MLCDWALAP')
    f4 = genes.Feature('CRISPR Cas1', (1210, 1400), '', 1, '', 4e-20, 'awesome protein', 'MCGYVERS')
    reannotated_operon = genes.Operon('operon1', '/path/to/contig.fa.gz', 0, 3400, [f1, f2, f3, f4])

    f1 = genes.Feature('CRISPR Cas1', (12, 400), '', 1, '', 4e-20, 'awesome protein', 'MCGYVERS')
    f2 = genes.Feature('RNA polymerase', (410, 601), '', 1, '', 4e-20, 'awesome protein', 'MRRQTHPR')
    f3 = genes.Feature('CRISPR Cas4 endonuclease', (620, 1205), '', 1, '', 4e-20, 'awesome protein', 'MLCDWALAP')
    f4 = genes.Feature('CRISPR Cas1', (1210, 1400), '', 1, '', 4e-20, 'awesome protein', 'MCGYVERS')
    reannotated_operon2 = genes.Operon('operon2', '/path/to/contig.fa.gz', 0, 3400, [f1, f2, f3, f4])

    actual = reannotation.count_cluster_reannotations([operon, operon2], {'operon1': reannotated_operon, 'operon2': reannotated_operon2})
    expected = {('cas1', 1): {'CRISPR Cas1': 2},
                ('cas2', 1): {'DNA polymerase': 1, 'RNA polymerase': 1},
                ('cas4', 1): {'CRISPR Cas4': 1, 'CRISPR Cas4 endonuclease': 1},
                ('cas1', 2): {'CRISPR Cas1': 2}}
    assert actual == expected
