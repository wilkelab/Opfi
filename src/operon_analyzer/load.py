import gzip
from typing import Optional, Iterator, IO
from Bio import SeqIO
from Bio.Seq import Seq
from operon_analyzer import genes, parse


def load_sequence(operon: genes.Operon) -> Optional[Seq]:
    """ Loads the DNA sequence for a given operon's contig. """
    with gzip.open(operon.contig_filename, 'rt') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            if record.id == operon.contig:
                return Seq(str(record.seq))


def load_gzipped_operons(filename: str) -> Iterator[genes.Operon]:
    with gzip.open(filename, 'rt') as f:
        for operon in load_operons(f):
            yield operon


def load_operons(handle: IO[str]) -> Iterator[genes.Operon]:
    yield from parse.assemble_operons(parse.read_pipeline_output(handle))
