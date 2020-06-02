import csv
from crisposon.output_writers import FIELDNAMES
from collections import defaultdict
from typing import Tuple, Iterator, IO
from operon_analyzer.genes import Feature, Operon


PipelineRecord = Tuple[str, str, str, str, str, str, str, str, str, str]
Coordinates = Tuple[int, int]


def assemble_operons(lines: Iterator[PipelineRecord]) -> Iterator[Operon]:
    """
    Takes the output from the CRISPR-transposon pipeline, loads all features,
    and assembles them into putative Operons.
    """
    next(lines)  # skip header
    operon_features = defaultdict(list)
    for line in lines:
        contig, coordinates, feature = _parse_feature(line)
        operon_features[(contig, coordinates)].append(feature)
    for (contig, (start, end)), features in operon_features.items():
        yield Operon(contig, start, end, features)


def parse_coordinates(coordinates: str) -> Coordinates:
    """ Parses base pair ranges denoting position in a contig. """
    start, end = coordinates.split("..")
    return int(start), int(end)


def read_pipeline_output(handle: IO[str]) -> Iterator[PipelineRecord]:
    """ Reads the CSV file produced by the CRISPR-transposon pipeline """
    reader = csv.reader(handle)
    for record in reader:
        if record == FIELDNAMES:
            continue
        yield record


def _parse_feature(line: PipelineRecord) -> Tuple[str, Coordinates, Feature]:
    """ Creates a Feature from a line of output from a CSVReader """
    contig = line[0]
    coordinates = parse_coordinates(line[1])
    feature = line[2]
    feature_coordinates = parse_coordinates(line[3])
    query_orfid = line[4]
    strand = int(line[5]) if line[5] else None
    hit_accession = line[6]
    hit_eval = float(line[7]) if line[7] else None
    description, sequence = line[8], line[9]
    return contig, coordinates, Feature(
        feature,
        feature_coordinates,
        query_orfid,
        strand,
        hit_accession,
        hit_eval,
        description,
        sequence)
