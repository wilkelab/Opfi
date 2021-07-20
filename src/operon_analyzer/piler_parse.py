""" Parses piler-cr output to extract the sequences of every spacer. """

from typing import List
from collections import namedtuple
from Bio.Seq import Seq


RepeatSpacer = namedtuple('RepeatSpacer', ['position', 'repeat_len', 'spacer_len', 'sequence'])
BrokenSpacer = namedtuple('BrokenSpacer', ['position', 'repeat_len', 'sequence'])


def parse_pilercr_output(text: str, start: int, end: int):
    """ Takes the text of :program:`pilercr` raw output, and a start and end coordinate for an operon's
    neighborhood, and extracts all spacers in that region. """
    if text is None:
        # There were no spacers in this contig
        return []
    text = text.split("\n")
    text = _find_first_entry(text)
    if text is None:
        # There were no spacers in this contig
        return []
    results = []
    while True:
        repeat_spacers, text = _parse_entry(text)
        valid_spacers = []
        for rs in repeat_spacers:
            if start <= rs.position <= end:
                valid_spacers.append(rs)
        if valid_spacers:
            results.append(valid_spacers)
        text = _find_next_entry(text)
        if not text:
            break
    return results


def _find_next_entry(text: List[str]):
    for n, line in enumerate(text):
        if line.startswith(">"):
            return text[n:]
        elif line.startswith("SUMMARY"):
            return False
    return False


def _find_first_entry(text: List[str]):
    for n, line in enumerate(text):
        if line.startswith(">"):
            return text[n:]


def _parse_entry(text: List[str]):
    repeat_spacers = []
    assert text[0].startswith(">")
    for n, line in enumerate(text[4:]):
        line = line.strip()
        if line.startswith("="):
            return repeat_spacers, text[n:]
        else:
            data = line.split()
            if len(data) == 7:
                position, repeat_len, _, spacer_len, _, _, spacer_seq = data
                if "N" in spacer_seq:
                    # don't allow spacers with any uncertainty
                    continue
                repeat_spacers.append(RepeatSpacer(int(position), int(repeat_len), int(spacer_len), Seq(spacer_seq)))
            elif len(data) == 6:
                position, repeat_len, _, _, _, broken_spacer_seq = data
                if "N" in broken_spacer_seq:
                    # don't allow spacers with any uncertainty
                    continue
                repeat_spacers.append(BrokenSpacer(int(position), int(repeat_len), Seq(broken_spacer_seq)))
