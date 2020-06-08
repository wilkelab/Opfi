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
    contig_filename = line[0]
    contig = line[1]
    coordinates = parse_coordinates(line[2])
    feature = line[3]
    feature_coordinates = parse_coordinates(line[4])
    query_orfid = line[5]
    strand = int(line[6]) if line[6] else None
    hit_accession = line[7]
    description = line[8]
    hit_eval = float(line[9]) if line[9] else None
    bit_score = float(line[10]) if line[10] else None

    # Additional blast alignment statistics
    raw_score = float(line[11]) if line[11] else None
    aln_len = int(line[12]) if line[12] else None
    pident = float(line[13]) if line[13] else None
    nident = int(line[14]) if line[14] else None
    mismatch = int(line[15]) if line[15] else None
    positive = int(line[16]) if line[16] else None
    gapopen = int(line[17]) if line[17] else None
    gaps = int(line[18]) if line[18] else None
    ppos = float(line[19]) if line[19] else None
    qcovhsp = int(line[20]) if line[20] else None

    # aligned part of query sequence
    sequence = line[21]

    if feature == "CRISPR array":
        copies, repeat, spacer = description.split(",")
        _, count = copies.split()
        feature = f"CRISPR array ({count})"
    return contig, coordinates, Feature(
        feature,
        feature_coordinates,
        query_orfid,
        strand,
        hit_accession,
        hit_eval,
        bit_score,
        description,
        sequence)
