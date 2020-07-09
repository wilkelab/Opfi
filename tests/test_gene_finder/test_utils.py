import pytest
from gene_finder.utils import get_neighborhood_ranges

def _build_hit_dictionary(coords):
    hits = {}
    for coord in coords:
        key = "hit_{}_{}".format(coord[0], coord[1])
        hits[key] = {}
        hits[key]["Query_start-pos"] = coord[0]
        hits[key]["Query_end-pos"] = coord[1]
    return hits


@pytest.mark.parametrize('hit_coords,expected_num_neighborhoods,expected_ranges', [
    ([(1,500), (3500, 3600), (6000, 6500)], 3, [(0, 1500), (2500, 4600), (5000, 7500)]),
    ([(1,500), (600, 11400), (20000, 20250)], 2, [(0, 12400), (19000, 21250)]),
    ([(1,500), (500, 1000), (1000, 1500)], 1, [(0, 2500)]),
    ([(500, 1), (600, 11400), (20000, 20250)], 2, [(0, 12400), (19000, 21250)]),
    ([(1,500), (500,1), (400, 600)], 1, [(0, 1600)]),
    ([(1, 500), (400, 2)], 1, [(0, 1500)]),
    ([(1, 500), (400, 1100), (1101, 1200)], 1, [(0, 2200)]),
    ([(2500, 2000)], 1, [(1000, 3500)])
    ])
def test_get_neighborhood_ranges(hit_coords, expected_num_neighborhoods, expected_ranges):
    hits = _build_hit_dictionary(hit_coords)
    neighborhoods = get_neighborhood_ranges(hits, 100000, span=1000)
    assert len(neighborhoods) == expected_num_neighborhoods
    for nbh, expected_range in zip(neighborhoods, expected_ranges):
        assert nbh == expected_range


@pytest.mark.parametrize('bait_ORF_coords,expected_last_range', [
    ([(6000, 6500)], (5000, 7500)),
    ([(1,500)], (0, 1500)),
    ([(98000, 99000)], (97000, 100000)),
    ([(98000, 99800)], (97000, 100000))
    ])
def test_get_neighborhood_ranges_bounds(bait_ORF_coords, expected_last_range):
    """
    Test that `get_neighborhood_ranges` always reports candidate region coordinates
    that are actually within the bounds of the parent contig.
    """
    hits = _build_hit_dictionary(bait_ORF_coords)
    neighborhoods = get_neighborhood_ranges(hits, contig_len=100000, span=1000)
    assert neighborhoods[0] == expected_last_range
