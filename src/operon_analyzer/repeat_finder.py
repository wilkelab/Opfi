# Identifies inverted and direct repeats and converts them to Feature objects
import subprocess
import tempfile
from operon_analyzer import genes, load
from Bio.Seq import Seq
from typing import Optional, Tuple, Iterator
from collections import namedtuple
import re


Repeat = namedtuple('Repeat', ['upstream_sequence', 'upstream_start', 'downstream_sequence', 'downstream_start'])
alignment_regex = re.compile(r'(?P<count>\d+)(?P<kind>\w)')


def find_inverted_repeats(operon: genes.Operon, buffer_around_operon: int, min_repeat_size: int):
    """
    Searches an operon and the DNA flanking it for inverted repeats. If found, they will be added to the operon
    as Feature objects with the name "TIR" and the sequence. The strand will be set to 1 for the upstream sequence and
    -1 for the downstream sequence.

    operon:     the Operon object
    buffer_around_operon:   the number of base pairs on either side of the operon to search in addition to the operon's internal sequence
    min_repeat_size:        the minimum number of base pairs that an inverted repeat must have
    """
    contig_sequence = load.load_sequence(operon)
    number_found = _find_inverted_repeats(operon, contig_sequence, buffer_around_operon, min_repeat_size)
    return number_found


def _find_inverted_repeats(operon: genes.Operon, contig_sequence: Seq, buffer_size: int, min_repeat_size: int):
    """ Performs the inverted repeat search on a given DNA sequence. """
    buffered_sequence = _get_buffered_operon_sequence(operon, contig_sequence, buffer_size)

    # We are looking for systems where the inverted repeats are outside of the operon we're interested in, so we set a spacer length that
    # precludes finding both inside of the bounds of the operon's features. It's still possible for one to be inside and one to be outside though.
    lower, upper = operon.feature_region
    min_spacer_size = upper - lower
    perfects, imperfects = _run_grf(buffered_sequence, min_spacer_size, min_repeat_size)

    # Go through each result, convert it to a pair of Features, and add the Features to the Operon
    for start, end, alignment in _parse_grf_results(perfects):
        upstream_feature, downstream_feature = _parse_repeats(start, end, alignment, buffered_sequence, buffer_size, operon.start, 'perfect')
        operon._features.append(upstream_feature)
        operon._features.append(downstream_feature)

    for start, end, alignment in _parse_grf_results(imperfects):
        corrected_start, corrected_end = _calculate_corrected_coordinates(start, end, buffer_size, operon.start)
        upstream_feature, downstream_feature = _parse_repeats(start, end, alignment, buffered_sequence, buffer_size, operon.start, 'imperfect')
        operon._features.append(upstream_feature)
        operon._features.append(downstream_feature)


def _calculate_corrected_coordinates(start, end, buffer_size, operon_start):
    """
    When setting the Feature coordinates, we need to account for all the adjustments we've made to the sequence.
    We also convert from GRF's 1-based indexes to 0-based.
    """
    corrected_start = start - buffer_size + operon_start - 1
    corrected_end = end - buffer_size + operon_start - 1
    return corrected_start, corrected_end


def _parse_repeats(start, end, alignment, buffered_sequence, buffer_size: int, operon_start: int, perfect_status: str):
    # Extract the raw sequences of the inverted repeats and convert them to an alignment string (i.e. add gaps for deletions/insertions)
    upstream_size, downstream_size = _parse_alignment_size(alignment)
    raw_upstream_seq = buffered_sequence[start - 1:start + upstream_size - 1]
    raw_downstream_seq = str(Seq(buffered_sequence[end - downstream_size: end]).reverse_complement())
    upstream_seq, downstream_seq = _format_aligned_sequences(raw_upstream_seq, raw_downstream_seq, alignment)

    # Figure out where the repeats are using the original coordinate system and create Feautre objects for each repeat
    corrected_start, corrected_end = _calculate_corrected_coordinates(start, end, buffer_size, operon_start)
    upstream_feature = _make_tir_feature(corrected_start, upstream_size, upstream_seq, 1, perfect_status)
    downstream_feature = _make_tir_feature(corrected_end - downstream_size + 1, downstream_size, downstream_seq, -1, perfect_status)
    return upstream_feature, downstream_feature


def _make_tir_feature(start: int, size: int, seq: str, strand: int, perfect_status: str) -> genes.Feature:
    """ Converts an inverted repeat into a Feature object. """
    assert len(seq) > 0
    assert start > 0
    assert size > 0
    assert strand in (-1, 1)
    return genes.Feature(f"TIR {seq}",
                         (start, start + size),
                         '',
                         strand,
                         '',
                         None,
                         f'{perfect_status}',
                         seq,
                         None)


def _get_buffered_operon_sequence(operon: genes.Operon, sequence: Seq, buffer_size: int) -> str:
    """ Selects a DNA sequence for an operon with a buffer on each side. """
    # operon.start and operon.end are 1-based indexes so we subtract 1 from operon.start
    start = max(0, operon.start - 1 - buffer_size)
    end = min(len(sequence), operon.end + buffer_size)
    sequence = str(sequence[start:end])
    return sequence


def _run_grf(sequence: str, min_spacer_size: int, min_repeat_size: int) -> Optional[Tuple[str, str]]:
    """
    Runs GenericRepeatFinder and returns the raw text of the perfect and imperfect spacers that are found.

    sequence: the sequence of the putative operon with buffers on each side
    min_repeat_size: the smallest inverted or direct repeat allowable. Must be >= 5
    """
    assert min_repeat_size >= 5, "min_repeat_size must be at least 5"
    assert len(sequence) > 0

    with tempfile.NamedTemporaryFile('w') as contig_file, tempfile.TemporaryDirectory() as grf_dir:
        contig_file.write(f">sequence\n{sequence}")
        # move the cursor to the beginning of the file or else we won't read any data from the temp file
        contig_file.seek(0)
        command = ["grf-main",
                   "-i", contig_file.name,
                   "-o", grf_dir,
                   "-c", "0",
                   "-f", "1",
                   "-s", str(min_repeat_size),
                   "--min_tr", str(min_repeat_size),
                   "--min_space", str(min_spacer_size),
                   "--max_space", str(len(sequence))]
        result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if result.returncode != 0:
            return None
        with open(f"{grf_dir}/perfect.spacer.id") as f:
            perfect = f.read().strip()
        with open(f"{grf_dir}/imperfect.id") as f:
            imperfect = f.read().strip()
        return perfect, imperfect


def _parse_grf_results(raw_text: str) -> Iterator[Tuple[int, int, str]]:
    """ Parses the raw text of GenericRepeatFinder output. Assumes that only
    IDs are written to the file, sequences should be omitted. If no results
    were found, the text will be empty. """
    if not raw_text:
        raise StopIteration
    for line in raw_text.split("\n"):
        yield _parse_repeat_id(line)


def _parse_repeat_id(repeat_id: str) -> Tuple[int, int, str]:
    """
    Parses one line from GenericRepeatFinder output.
    start is the location of the beginning of the upstream repeat.
    end is the location of the end of the downstream repeat.
    alignment is a serialization of the difference between the two sequences
    (see _parse_alignment_size for details).
    """
    _, start, end, alignment = repeat_id.split(":")
    return int(start), int(end), alignment


def _parse_alignment_size(alignment: str) -> Tuple[int, int]:
    """
    We are given the difference (if any) between the upstream and downstream
    repeat in a serialized code. In order to determine the actual sequence of the repeat,
    we have to parse this string just to determine the size of the repeats.

    The letters mean:
    m:  perfect match
    M:  aligned mismatch
    I:  insertion in upstream repeat
    D:  deletion in upstream repeat

    For example, `5m3I16m` would be a repeat with 5 homologous base pairs, 3 insertions
    in the upstream repeat relative to the downstream repeat, followed by 16 homologous
    base pairs.
    """
    sizes = {'m': (1, 1), 'M': (1, 1), 'I': (1, 0), 'D': (0, 1)}
    matches = alignment_regex.findall(alignment)
    upstream_total, downstream_total = 0, 0
    for count, kind in matches:
        count = int(count)
        upfactor, downfactor = sizes[kind]
        upstream_total += upfactor * count
        downstream_total += downfactor * count
    return upstream_total, downstream_total


def _format_aligned_sequences(raw_upstream_sequence: str,
                              raw_downstream_sequence: str,
                              alignment: str) -> Tuple[str, str]:
    """ Takes the raw sequences of each repeat and inserts dashes to indicate gaps, as appropriate. """
    matches = alignment_regex.findall(alignment)
    upseq = []
    downseq = []
    up_position = 0
    down_position = 0
    for count, kind in matches:
        count = int(count)
        if kind in ('m', 'M'):
            # mismatches don't require a dash
            upseq.append(raw_upstream_sequence[up_position: up_position + count])
            downseq.append(raw_downstream_sequence[down_position: down_position + count])
            up_position += count
            down_position += count
        elif kind == 'I':
            # insertions mean the upstream repeat has an extra base relative to the downstream repeat
            upseq.append(raw_upstream_sequence[up_position: up_position + count])
            downseq.append('-' * count)
            up_position += count
        elif kind == 'D':
            # deletions mean the downstream repeat has an extra base relative to the upstream repeat
            upseq.append('-' * count)
            downseq.append(raw_downstream_sequence[down_position: down_position + count])
            down_position += count
        else:
            raise ValueError(f"unexpected alignment value: {alignment}")
    # convert the lists that are accumulating strings to a single aligned sequence
    return "".join(upseq), "".join(downseq)
