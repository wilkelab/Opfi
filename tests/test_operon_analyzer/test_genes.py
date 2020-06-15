from operon_analyzer.genes import Feature, Operon


def test_feature_length():
    f = Feature('gene', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    assert len(f) == 100


def test_get_unique():
    f1 = Feature('gene1', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f2 = Feature('gene2', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f3 = Feature('gene3', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f4 = Feature('gene3', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    op = Operon('asdf', 0, 1000, [f1, f2, f3, f4])
    funique = op.get_unique('gene1')
    assert funique is not None


def test_get_unique_not_unique():
    f1 = Feature('gene1', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f2 = Feature('gene2', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f3 = Feature('gene3', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    f4 = Feature('gene3', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    op = Operon('asdf', 0, 1000, [f1, f2, f3, f4])
    fnotunique = op.get_unique('gene3')
    assert fnotunique is None


def test_feature_repr():
    f = Feature('gene1', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    label = str(f)
    assert label == '<Feature gene1 1..100>'
