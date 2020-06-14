from operon_analyzer.genes import Feature, Operon


def test_feature_length():
    f = Feature('gene', (1, 100), '', 1, '', 4e-12, '', 'MGWRN', 1234.9)
    assert len(f) == 100
