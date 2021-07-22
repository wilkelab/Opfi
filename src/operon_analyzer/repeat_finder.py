# Identifies inverted repeats and converts them to Feature objects
import re
import subprocess
import tempfile
from collections import namedtuple
from typing import Iterator, Optional, Tuple

from Bio.Seq import Seq

from operon_analyzer import genes, load

Repeat = namedtuple('Repeat', ['upstream_sequence', 'upstream_start', 'downstream_sequence', 'downstream_start'])
GRFResult = namedtuple('GRFResult', ['start', 'end', 'alignment'])
alignment_regex = re.compile(r'(?P<count>\d+)(?P<kind>\w)')


class BufferedSequence(object):
    """ Provides access to an operon's sequence, with a buffer flanking each side. This allows
    us to search for inverted repeats adjacent to a putative operon and get their exact sequences.
    All coordinates are internally converted to Python indicies so that we can easily slice sequences correctly,
    but we report the 1-based indexes to stay compatible with the rest of the library.
    """

    def __init__(self, operon: genes.Operon, contig_sequence: Seq, buffer_size: int):
        self._operon_start, self._operon_end = operon.feature_region
        self._contig_sequence = contig_sequence
        self._buffer_size = buffer_size

    @property
    def start(self) -> int:
        """ Returns the 0-based index of the beginning of the buffered sequence. """
        return max(0, self._operon_start - self._buffer_size)

    @property
    def end(self) -> int:
        """ Returns the 0-based, exclusive index of the end of the buffered sequence. """
        return min(len(self._contig_sequence), self._operon_end + self._buffer_size)

    @property
    def sequence(self) -> str:
        return self._contig_sequence[self.start:self.end]

    def __getitem__(self, key) -> str:
        """ Get slices of the contig sequence using coordinates that define the 0 position as the
        beginning of the buffered sequence, NOT the contig sequence. Used to get the sequence of
        the inverted repeats. """
        return self.sequence[key.start: key.stop]

    @property
    def operon_length(self) -> int:
        """ Gives the length of the putative operon. """
        return self._operon_end - self._operon_start + 1


def find_inverted_repeats(operon: genes.Operon, buffer_around_operon: int, min_repeat_size: int):
    """
    Searches an operon and the DNA flanking it for inverted repeats using :program:`GenericRepeatFinder`. 
    If found, they will be added to the operon as Feature objects with the name "TIR" and the sequence. 
    The strand will be set to 1 for the upstream sequence and -1 for the downstream sequence.

    Args:
        operon (Operon): The :class:`operon_analyzer.genes.Operon` object.
        buffer_around_operon (int): The number of base pairs on either side of the operon to search in addition to the operon's 
            internal sequence.
        min_repeat_size (int): The minimum number of base pairs that an inverted repeat must have.
    """
    contig_sequence = load.load_sequence(operon)
    _find_inverted_repeats(operon, contig_sequence, buffer_around_operon, min_repeat_size)


def _find_inverted_repeats(operon: genes.Operon, contig_sequence: Seq, buffer_size: int, min_repeat_size: int):
    """ Performs the inverted repeat search on a given DNA sequence. """
    buffered_sequence = BufferedSequence(operon, contig_sequence, buffer_size)

    # We are looking for systems where the inverted repeats are outside of the operon we're interested in, so we set a
    # spacer length that precludes finding both inside of the bounds of the operon's features. It's still possible for
    # one to be inside and one to be outside though.
    perfects, imperfects = _run_grf(buffered_sequence, min_repeat_size)

    # Go through each result, convert it to a pair of Features, and add the Features to the Operon
    n = 0  # must set in case we don't have any perfect results
    for n, result in enumerate(_parse_grf_results(perfects)):
        upstream_feature, downstream_feature = _parse_repeats(result, buffered_sequence, n)
        operon._features.append(upstream_feature)
        operon._features.append(downstream_feature)

    count = n
    for n, result in enumerate(_parse_grf_results(imperfects)):
        upstream_feature, downstream_feature = _parse_repeats(result, buffered_sequence, n + count + 1)
        operon._features.append(upstream_feature)
        operon._features.append(downstream_feature)


def _parse_repeats(result: GRFResult, buffered_sequence: BufferedSequence, number: int):
    # Extract the raw sequences of the inverted repeats and convert them to an alignment string (i.e. add gaps for deletions/insertions)
    upstream_size, downstream_size = _parse_alignment_size(result.alignment)
    raw_upstream_seq = str(buffered_sequence[result.start: result.start + upstream_size])
    raw_downstream_seq = str(buffered_sequence[result.end - downstream_size: result.end].reverse_complement())
    upstream_alignment, downstream_alignment = _format_aligned_sequences(raw_upstream_seq, raw_downstream_seq, result.alignment)

    # Figure out where the repeats are using the original coordinate system and create Feautre objects for each repeat
    upstream_feature = _make_tir_feature(buffered_sequence.start + result.start, upstream_size, raw_upstream_seq, upstream_alignment, 1, number)
    downstream_feature = _make_tir_feature(buffered_sequence.start + result.end - downstream_size, downstream_size, raw_downstream_seq, downstream_alignment, -1, number)
    return upstream_feature, downstream_feature


def _make_tir_feature(start: int, size: int, seq: str, alignment: str, strand: int, number) -> genes.Feature:
    """ Converts an inverted repeat into a Feature object. """
    assert len(seq) > 0
    assert start >= 0
    assert size > 0
    assert strand in (-1, 1)
    return genes.Feature(f"IR #{number} ({len(alignment)} bp)",
                         (start, start + size),
                         '',
                         strand,
                         '',
                         None,
                         f'{seq} {alignment}',
                         seq,
                         None)


def _run_grf(buffered_sequence: BufferedSequence, min_repeat_size: int) -> Optional[Tuple[str, str]]:
    """
    Runs GenericRepeatFinder and returns the raw text of the perfect and imperfect spacers that are found.

    sequence: the sequence of the putative operon with buffers on each side
    min_repeat_size: the smallest inverted or direct repeat allowable. Must be >= 5
    """
    assert min_repeat_size >= 5, "min_repeat_size must be at least 5"

    with tempfile.NamedTemporaryFile('w') as contig_file, tempfile.TemporaryDirectory() as grf_dir:
        contig_file.write(f">sequence\n{buffered_sequence.sequence}")
        # move the cursor to the beginning of the file or else we won't read any data from the temp file
        contig_file.seek(0)
        command = ["grf-main",
                   "-i", contig_file.name,
                   "-o", grf_dir,
                   "-c", "0",
                   "-f", "1",
                   "-s", str(min_repeat_size),
                   "--min_tr", str(min_repeat_size),
                   "--min_space", str(buffered_sequence.operon_length),
                   "--max_space", str(len(buffered_sequence.sequence))]
        result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if result.returncode != 0:
            return None
        with open(f"{grf_dir}/perfect.spacer.id") as f:
            perfect = f.read().strip()
        with open(f"{grf_dir}/imperfect.id") as f:
            imperfect = f.read().strip()
        return perfect, imperfect


def _parse_grf_results(raw_text: str) -> Iterator[GRFResult]:
    """ Parses the raw text of GenericRepeatFinder output. Assumes that only
    IDs are written to the file, sequences should be omitted. If no results
    were found, the text will be empty. """
    if not raw_text:
        return None
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
    # GRF uses 1-based indexes, with inclusive ends. Subtracting 1 from the start
    # coordinate allows us to use the coordinates as slice indices
    result = GRFResult(int(start) - 1, int(end), alignment)
    assert result.start >= 0
    return result


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
