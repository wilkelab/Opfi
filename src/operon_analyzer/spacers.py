import gzip
from typing import List, Union, Optional, Tuple
from collections import namedtuple
import tempfile
import subprocess
import statistics
import json
import multiprocessing
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from operon_analyzer import genes
import piler_parse
import parasail


MATCH_SCORE = 5
MISMATCH_PENALTY = -1
GAP_OPEN_PENALTY = -8
GAP_EXTEND_PENALTY = -8

AlignmentResult = namedtuple('AlignmentResult', ['match_count', 'spacer_sequence', 'contig_sequence', 'spacer_alignment', 'contig_alignment', 'contig_start', 'contig_end', 'strand'])
MAX_SPACER_LENGTH_BP = 40  # Assume spacers longer than this were incorrectly parsed and should be ignored
Spacer = Union[piler_parse.RepeatSpacer, piler_parse.BrokenSpacer]
Array = List[Spacer]


def find_self_targeting_spacers(operon: genes.Operon, min_matches: int, num_processes: int = 32) -> List[genes.Feature]:
    """
    Finds all self-targeting spacers in an Operon's CRISPR arrays,
    and returns Feature objects for each of those (with the name "CRISPR target")

    Steps:
      - Loads the sequence of the contig that the operon came from
      - Runs pilercr again and extracts the spacer sequences
      - Performs a pairwise local alignment of each spacer and the contig
      - For spacers above the match threshold, creates a Feature object
      - Returns the Features
    """
    contig_sequence = load_sequence(operon)
    tempseq = str(contig_sequence)
    tempseq = tempseq[:30] + 'ATTGAGCAACACCGAAATATTAACTATTAATTTTTCAACAAGCAGACAT' + tempseq[30:]
    assert contig_sequence is not None, f"Operon's sequence file cannot be found: {operon.contig_filename}."
    piler_data = run_piler(contig_sequence)
    arrays = piler_parse.parse_pilercr_output(piler_data)
    fixed_arrays = fix_arrays(arrays, contig_sequence)
    return run_alignments_parallel(fixed_arrays, contig_sequence, num_processes, min_matches)


def run_alignments_parallel(arrays: List[Array], contig: str, num_processes: int, min_matches: int):
    spacers = [spacer for array in arrays for spacer in array]
    pool = multiprocessing.Pool(min(num_processes, len(spacers)))
    results = []
    for spacer in spacers:
        result = pool.apply_async(align_spacer_to_contig, args=(spacer, contig, min_matches))
        results.append(result)
    pool.close()
    pool.join()

    features = []
    for result in results:
        ar = result.get()
        if ar:
            feature = build_feature_from_alignment(ar)
            features.append(feature)
    return features


def build_feature_from_alignment(ar: AlignmentResult) -> genes.Feature:
    description = json.dumps({"match_count": ar.match_count,
                              "spacer_alignment": ar.spacer_alignment,
                              "contig_alignment": ar.contig_alignment,
                              "spacer_sequence": str(ar.spacer_sequence),
                              "contig_sequence": str(ar.contig_sequence)})
    feature = genes.Feature('CRISPR target',
                            (ar.contig_start, ar.contig_end),
                            '',
                            ar.strand,
                            '',
                            None,
                            description,
                            '',
                            None)
    return feature


def align_spacer_to_contig(spacer: Seq,
                           contig: str,
                           min_matches: int):
    alignment_result = perform_local_pairwise_alignment(spacer, contig)
    if alignment_result.match_count < min_matches:
        return None
    return alignment_result


def find_target(spacer, contig):
    result = parasail.ssw(spacer, contig, -GAP_EXTEND_PENALTY, -GAP_OPEN_PENALTY, parasail.blosum62)
    target = contig[result.ref_begin1:result.ref_end1 + 1]
    top, bottom = pairwise2.align.localms(target, spacer, MATCH_SCORE, MISMATCH_PENALTY, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY)[0][:2]
    matches = sum([1 for t, b in zip(top, bottom) if t == b and t != '-'])
    return top, bottom, matches, result.ref_begin1, result.ref_end1, target


def align_contig_and_spacer(contig: Seq, spacer: Seq, spacer_location_start: int, spacer_location_end: int):
    # block out the actual spacer since it's by definition a perfect match
    censored_contig = str(contig[:spacer_location_start]) + "N"*len(spacer) + str(contig[spacer_location_end:])
    contig_alignment, spacer_alignment, match_count, contig_start, contig_end, target = find_target(spacer, censored_contig)
    return match_count, contig_alignment, spacer_alignment, contig_start, contig_end, target


def perform_local_pairwise_alignment(spacer: piler_parse.RepeatSpacer, contig: Seq) -> Tuple[int, int, str, str, int, int]:
    best_match_count = 0
    best_strand = None
    spacer_seqs = spacer.sequence, spacer.sequence.reverse_complement()

    for n, spacer_seq in enumerate(spacer_seqs):
        spacer_start = spacer.position + spacer.repeat_len
        spacer_end = spacer.position + spacer.repeat_len + spacer.spacer_len + 1
        results = align_contig_and_spacer(contig, str(spacer_seq), spacer_start, spacer_end)
        match_count = results[0]
        if match_count > best_match_count:
            best_match_count = match_count
            best_results = results
            best_strand = n

    match_count, contig_alignment, spacer_alignment, contig_start, contig_end, target = best_results
    strand = [1, -1][best_strand]

    return AlignmentResult(match_count,
                           str(spacer_seqs[best_strand]),
                           target,
                           spacer_alignment,
                           contig_alignment,
                           contig_start,
                           contig_end,
                           strand)


def run_piler(sequence: str) -> str:
    """ Runs pilercr on the given Seq sequence and returns the raw file output """
    with tempfile.NamedTemporaryFile('w') as contig_file, tempfile.NamedTemporaryFile('w') as pilercr_file:
        contig_file.write(f">sequence\n{sequence}")
        command = ["pilercr", "-in", contig_file.name,
                              "-out", pilercr_file.name,
                              "-minarray", "2",
                              "-quiet",
                              "-noinfo"]
        subprocess.call(command)
        with open(pilercr_file.name) as f:
            return f.read()


def fix_arrays(arrays: List[Array], contig: Seq) -> List[Array]:
    fixed_arrays = []
    for array in arrays:
        median_length = find_median_spacer_length(array)
        if not median_length:
            continue
        fixed_array = fix_broken_spacers(array, median_length, contig)
        fixed_arrays.append(fixed_array)
    return fixed_arrays


def find_median_spacer_length(array: List[Spacer]) -> Optional[int]:
    good_spacers = [rs.spacer_len for rs in array if type(rs) == piler_parse.RepeatSpacer]
    if not good_spacers:
        return None
    return int(statistics.median(good_spacers))


def fix_broken_spacers(array: List[Spacer], median_length: int, sequence: Seq) -> List[piler_parse.RepeatSpacer]:
    """ Pilercr sometimes reports that the last spacer in an array is much shorter than the rest.
    We assume that it should at least have the median length of the rest of the spacers and extend
    it to that length.
    """
    fixed_array = []
    for rs in array:
        # Figure out what the spacer sequence is
        if type(rs) is piler_parse.BrokenSpacer:
            spacer = str(sequence)[rs.position + rs.repeat_len - 1:rs.position + rs.repeat_len + median_length - 1]
            if len(spacer) <= MAX_SPACER_LENGTH_BP:
                fixed_array.append(piler_parse.RepeatSpacer(rs.position, rs.repeat_len, len(spacer), Seq(spacer)))
        else:
            fixed_array.append(rs)
    return fixed_array


def load_sequence(operon: genes.Operon):
    with gzip.open(operon.contig_filename, 'rt') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            if record.id == operon.contig:
                return Seq(str(record.seq))
