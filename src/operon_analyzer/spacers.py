from typing import List, Union, Optional, Tuple
from collections import namedtuple
import tempfile
import subprocess
import statistics
import json
import multiprocessing
from Bio.Seq import Seq
from operon_analyzer import genes, piler_parse, load
import parasail
import re


GAP_OPEN_PENALTY = 8
GAP_EXTEND_PENALTY = 8

AlignmentResult = namedtuple('AlignmentResult', ['match_count', 'spacer_sequence', 'contig_sequence', 'contig_start', 'contig_end', 'spacer_alignment', 'contig_alignment', 'comp_string', 'strand', 'spacer_order', 'array_length'])
MAX_SPACER_LENGTH_BP = 40  # Assume spacers longer than this were incorrectly parsed and should be ignored
Spacer = Union[piler_parse.RepeatSpacer, piler_parse.BrokenSpacer]
Array = List[Spacer]

cigar_regex = re.compile(r"(\d+)=")
start_regex = re.compile(r'^(\d+)D')


def find_self_targeting_spacers(operons: List[genes.Operon], min_matching_fraction: float, num_processes: int = 32):
    """
    For each given Operon, this will determine if CRISPR spacers target a location in the operon's parent contig.
    Matches with more than ``min_matching_fraction`` homology will be added to the Operon as a Feature named "CRISPR target".
    """
    pool = multiprocessing.Pool(min(num_processes, len(operons)))
    results = []
    for operon in operons:
        result = pool.apply_async(_align_operon_spacers, args=(operon, min_matching_fraction))
        results.append(result)
    pool.close()
    pool.join()
    for result in results:
        yield result.get()


def _align_operon_spacers(operon: genes.Operon, min_matching_fraction: float):
    """
    Aligns every spacer in an Operon to its parent contig. Any spacers with a target that has
    at least min_matching_fraction homology will be converted to Features with the name
    "CRISPR target" and added to the Operon object.
    """
    contig_sequence = load.load_sequence(operon)
    assert contig_sequence is not None, f"Operon's sequence file cannot be found: {operon.contig_filename}."
    spacers = _get_operon_spacers(operon.start, operon.end, contig_sequence)
    if not spacers:
        return operon

    # Look up where the arrays are. If a spacer targets the array it comes from, we don't create a CRISPR target.
    # It's either an assembly artifact or a duplication of the array.
    arrays = [feature for feature in operon.all_features if feature.name == 'CRISPR array']

    for spacer, spacer_order, array_length in spacers:
        ar = _perform_local_pairwise_alignment(spacer, spacer_order, array_length, str(contig_sequence))
        if not ar:
            continue
        matching_fraction = ar.match_count / len(ar.spacer_sequence)
        if matching_fraction >= min_matching_fraction:
            feature = _build_feature_from_alignment(ar, spacer)

            # Make sure the spacer doesn't target its own array
            for array in arrays:
                if array.start <= feature.start <= array.end and array.start <= spacer.position <= array.end:
                    # The spacer and target are both in the same array
                    break
            else:
                # The target is somewhere outside of the array where its spacer is
                operon._features.append(feature)
    return operon


def _get_operon_spacers(operon_start: int, operon_end: int, contig_sequence: Seq):
    """
    Runs pilercr and extracts all the spacer sequences found within the operon.
    Also adds the position of the spacer in the array, and the total array length.
    """
    piler_data = _run_piler(contig_sequence)
    arrays = piler_parse.parse_pilercr_output(piler_data, operon_start, operon_end)
    fixed_arrays = _fix_arrays(arrays, contig_sequence)
    spacers = []
    for array in fixed_arrays:
        for n, spacer in enumerate(array):
            if len(spacer.sequence) > MAX_SPACER_LENGTH_BP:
                continue
            spacers.append((spacer, n, len(array)))
    return spacers


def _perform_local_pairwise_alignment(spacer: piler_parse.RepeatSpacer,
                                      spacer_order: int,
                                      array_length: int,
                                      contig: Seq) -> Optional[AlignmentResult]:
    """
    Aligns a spacer and its reverse complement to the given contig and returns the better of the two matches.
    The spacer sequence itself is removed from the contig so that we don't just return the spacer itself.
    """
    best_score = 0
    best_strand = None
    assert str(spacer.sequence) in contig
    censored_contig = _build_censored_contig(spacer, contig)

    spacer_seqs = spacer.sequence, spacer.sequence.reverse_complement()

    for n, spacer_seq in enumerate(spacer_seqs):
        results = _align(str(spacer_seq), censored_contig)
        score = results[0]
        if score > best_score:
            best_score = score
            best_results = results
            best_strand = n

    _, match_count, contig_target, contig_start, contig_end, contig_alignment, spacer_alignment, comp_string = best_results
    strand = [1, -1][best_strand]

    alignment_result = AlignmentResult(match_count,
                                       str(spacer_seqs[best_strand]),
                                       contig_target,
                                       contig_start,
                                       contig_end,
                                       spacer_alignment,
                                       contig_alignment,
                                       comp_string,
                                       strand,
                                       spacer_order,
                                       array_length)

    # We shouldn't have more matching base pairs than base pairs
    assert alignment_result.match_count <= len(alignment_result.spacer_sequence), str(alignment_result)
    assert alignment_result.match_count <= len(alignment_result.contig_sequence), str(alignment_result)

    if "N" in alignment_result.contig_sequence:
        # don't permit contigs with any uncertainty
        return None
    return alignment_result


def _build_censored_contig(spacer: piler_parse.RepeatSpacer, contig: Seq) -> str:
    """
    Replace the spacer sequence itself in the contig with repeating "N"s so that
    the local alignment doesn't just find the spacer itself.
    """
    spacer_start = spacer.position + spacer.repeat_len
    spacer_end = spacer.position + spacer.repeat_len + spacer.spacer_len
    censored = str(contig[:spacer_start]) + "N"*spacer.spacer_len + str(contig[spacer_end:])
    assert len(censored) == len(contig), f"{len(contig)} {len(censored)}"
    return censored


def _align(spacer: str, contig: str) -> Tuple[float, int, str, int, int, str, str, str]:
    """ Run the local pairwise alignment of two strings and return alignment data. """
    # perform the alignment
    result = parasail.sw_trace(spacer, contig, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, parasail.nuc44)

    # extract pertinent data from the alignment result
    cigar_text = result.cigar.decode.decode("utf8")
    match_count = _count_cigar_matches(cigar_text)
    contig_target_sequence = result.traceback.ref.replace("-", "")
    contig_end = result.end_ref
    contig_start = contig_end - len(contig_target_sequence) + 1
    return result.score, match_count, contig_target_sequence, contig_start, contig_end, result.traceback.ref, result.traceback.query, result.traceback.comp


def _build_feature_from_alignment(ar: AlignmentResult, spacer: piler_parse.RepeatSpacer) -> genes.Feature:
    """ Alignment results are converted to Features so that we can plot them easily. """
    description = json.dumps({"match_count": ar.match_count,
                              "matching_fraction": (ar.match_count / len(str(ar.spacer_sequence))),
                              "spacer_alignment": ar.spacer_alignment,
                              "contig_alignment": ar.contig_alignment,
                              "contig_start": ar.contig_start,
                              "contig_end": ar.contig_end,
                              "spacer_sequence": str(ar.spacer_sequence),
                              "contig_sequence": str(ar.contig_sequence),
                              "comp_string": ar.comp_string,
                              "spacer_order": ar.spacer_order,
                              "spacer_position": spacer.position,
                              "array_length": ar.array_length})
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


def _count_cigar_matches(string: str) -> int:
    """ parasail provides a CIGAR string to encode the alignment. We parse this to determine the number of exact matches. """
    matches = cigar_regex.findall(string)
    return sum([int(match) for match in matches])


def _run_piler(sequence: str) -> str:
    """ Runs pilercr on the given Seq sequence and returns the raw file output """
    with tempfile.NamedTemporaryFile('w', dir='/tmp') as contig_file, tempfile.NamedTemporaryFile('w', dir='/tmp') as pilercr_file:
        contig_file.write(f">sequence\n{sequence}")
        command = ["pilercr", "-in", contig_file.name,
                              "-out", pilercr_file.name,
                              "-minarray", "2",
                              "-quiet",
                              "-noinfo"]
        # Piler crashes on some inputs. We pipe stderr to suppress the error messages.
        # We may be losing results when this occurs, but it's pretty rare and there's
        # nothing we can do about it either way.
        result = subprocess.call(command, stderr=subprocess.DEVNULL)
        if result != 0:
            return None
        with open(pilercr_file.name) as f:
            return f.read()


def _fix_arrays(arrays: List[Array], contig: Seq) -> List[Array]:
    """
    Ensures that the final spacer in each array is a reasonable length. For some reason, they
    are occasionally substantially shorter than all preceeding spacers in the array.
    """
    fixed_arrays = []
    for array in arrays:
        median_length = _find_median_spacer_length(array)
        if not median_length:
            continue
        fixed_array = _fix_array(array, median_length, contig)
        fixed_arrays.append(fixed_array)
    return fixed_arrays


def _find_median_spacer_length(array: List[Spacer]) -> Optional[int]:
    good_spacers = [rs.spacer_len for rs in array if type(rs) == piler_parse.RepeatSpacer]
    if not good_spacers:
        return None
    return int(statistics.median(good_spacers))


def _fix_array(array: List[Spacer], median_length: int, sequence: Seq) -> List[piler_parse.RepeatSpacer]:
    """
    Pilercr sometimes reports that the last spacer in an array is much shorter than the rest.
    We assume that it should at least have the median length of the rest of the spacers and extend
    it to that length.
    """
    fixed_array = []
    for rs in array:
        # Figure out what the spacer sequence is
        if type(rs) is piler_parse.BrokenSpacer:
            fixed_rs = _fix_broken_spacer(rs, median_length, sequence)
            if fixed_rs:
                fixed_array.append(fixed_rs)
        else:
            fixed_array.append(rs)
    return fixed_array


def _fix_broken_spacer(rs: piler_parse.BrokenSpacer, median_length: int, sequence: Seq) -> Optional[piler_parse.RepeatSpacer]:
    spacer = str(sequence)[rs.position + rs.repeat_len: rs.position + rs.repeat_len + median_length]
    if len(spacer) > MAX_SPACER_LENGTH_BP:
        return None
    return piler_parse.RepeatSpacer(rs.position, rs.repeat_len, len(spacer), Seq(spacer))
