# Identifies inverted and direct repeats and converts them to Feature objects
import subprocess
import tempfile
from typing import Optional


def _run_grf(operon_sequence: str, buffer_size: int, min_repeat_size=5, num_threads=4) -> Optional[str]:
    """
    Runs GenericRepeatFinder and returns the raw text of the perfect and imperfect spacers that are found.

    operon_sequence: the sequence of the putative operon with buffers on each side
    buffer_size:     the amount of extra DNA on each side of the putative operon. If unequal, this should be the larger of the two.
    min_repeat_size: the smallest inverted or direct repeat allowable. Must be >= 5
    num_threads:     number of threads to use
    """
    assert min_repeat_size >= 5, "min_repeat_size must be at least 5"
    assert len(operon_sequence) > 0
    with tempfile.NamedTemporaryFile('w', dir='/tmp') as contig_file, tempfile.TemporaryDirectory(dir='/tmp') as grf_dir:
        contig_file.write(f">sequence\n{operon_sequence}")
        command = ["grf-main", "-i", contig_file.name,
                   "-o", grf_dir,
                   "-c", "0",
                   "--min_tr", str(min_repeat_size),
                   "-s", str(min_repeat_size),
                   "-t", str(num_threads),
                   "--min_space", str(len(operon_sequence) - buffer_size),
                   "--max_space", str(len(operon_sequence))]

        # Piler crashes on some inputs. We pipe stderr to suppress the error messages.
        # We may be losing results when this occurs, but it's pretty rare and there's
        # nothing we can do about it either way.
        result = subprocess.call(command, stderr=subprocess.DEVNULL)
        if result != 0:
            return None
        with open(f"{grf_dir}/perfect.spacer.fasta") as f:
            perfect = f.read()
        with open(f"{grf_dir}/imperfect.fasta") as f:
            imperfect = f.read()
        return perfect, imperfect
