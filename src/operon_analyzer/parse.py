import csv
from gene_finder.output_writers import FIELDNAMES
from typing import Tuple, Iterator, IO
from operon_analyzer.genes import Feature, Operon
import sys


csv.field_size_limit(sys.maxsize)

PipelineRecord = Tuple[str, str, str, str, str, str, str, str, str, str]
Coordinates = Tuple[int, int]


def load_operons(handle: IO[str]) -> Iterator[Operon]:
    yield from assemble_operons(read_pipeline_output(handle))


def assemble_operons(lines: Iterator[PipelineRecord]) -> Iterator[Operon]:
    """
    Takes the output from the CRISPR-transposon pipeline, loads all features,
    and assembles them into putative Operons.

    To keep things memory efficient while not allowing redundant operons from being loaded,
    we keep track of the hash of each operon, which is just an integer. Even for very
    large metagenomic databases this should use less than a gigabyte of memory.
    """
    hashes = set()
    first = next(lines)
    contig, contig_filename, coordinates, feature = _parse_feature(first)
    current_contig = (contig, contig_filename, coordinates)
    features = [feature]

    for line in lines:
        contig, contig_filename, coordinates, feature = _parse_feature(line)
        if (contig, contig_filename, coordinates) != current_contig:
            (current_contig, current_contig_filename, (current_start, current_end)) = current_contig
            operon = Operon(current_contig, current_contig_filename, current_start, current_end, features)
            h = hash(operon)
            if h not in hashes:
                hashes.add(h)
                yield operon
            features = [feature]
            current_contig = (contig, contig_filename, coordinates)
        else:
            features.append(feature)
    (current_contig, current_contig_filename, (current_start, current_end)) = current_contig
    operon = Operon(current_contig, current_contig_filename, current_start, current_end, features)
    h = hash(operon)
    if h not in hashes:
        hashes.add(h)
        yield operon


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

    # Piler-cr and BLAST both use 1-based indices, but Opfi uses 0-based indices.
    # To make both coordinate systems consistent, we subtract 1 from the start
    # since feature coordinates come directly from those tools.
    # If features are on the reverse strand, the second coordinate will be larger
    # than the first, but operon_analyzer assumes the start is always less than the 
    # end
    first_coord, second_coord = parse_coordinates(line[3])
    feature_start = min(first_coord, second_coord) - 1
    feature_end = max(first_coord, second_coord)

    query_orfid = line[4]
    strand = int(line[5]) if line[5] else None
    hit_accession = line[6]
    hit_eval = float(line[7]) if line[7] else None
    description = line[8]
    sequence = line[9]
    bit_score = None
    # A newer version of the pipeline writes output with
    # an additonal 12 columns.
    # We'll grab these values too if they're present, but only
    # the bitscore (column 10) is actually used for now
    if len(line) == 22:
        # additional blast alignment stats
        bit_score = float(line[10]) if line[10] else None
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

        # fasta file name
        contig_filename = line[21] if line[21] else ''
    else:
        contig_filename = ''

    return contig, contig_filename, coordinates, Feature(
        feature,
        (feature_start, feature_end),
        query_orfid,
        strand,
        hit_accession,
        hit_eval,
        description,
        sequence,
        bit_score)
