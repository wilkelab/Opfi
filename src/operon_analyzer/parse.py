from collections import defaultdict
from typing import Tuple, List, Iterator
from operon_analyzer.genes import Feature, Operon


def parse_pipeline_results(lines: Iterator[Tuple]) -> Iterator[Operon]:
    """
    Takes the output from the CRISPR-transposon pipeline, loads all features,
    and assembles them into putative Operons.
    """
    operon_features = defaultdict(list)
    for line in lines:
        contig, coordinates, feature = _parse_feature(line)
        operon_features[(contig, coordinates)].append(feature)
    for (contig, (start, end)), features in operon_features.items():
        yield Operon(contig, start, end, features)


def _parse_coordinates(coordinates: str) -> Tuple[int, int]:
    """ Parses base pair ranges denoting position in a contig. """
    start, end = coordinates.split("..")
    return int(start), int(end)


def _parse_feature(line: List[str]) -> (str, Tuple[int, int], Feature):
    """ Creates a Feature from a line of output from a CSVReader """
    contig = line[0]
    coordinates = _parse_coordinates(line[1])
    feature = line[2]
    feature_coordinates = _parse_coordinates(line[3])
    query_orfid = line[4]
    strand = int(line[5]) if line[5] else None
    hit_accession = line[6]
    hit_eval = float(line[7]) if line[7] else None
    description, sequence = line[8:]
    return contig, coordinates, Feature(
        feature,
        feature_coordinates,
        query_orfid,
        strand,
        hit_accession,
        hit_eval,
        description,
        sequence)
