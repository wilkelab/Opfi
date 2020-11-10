import gzip
from typing import Optional
from Bio import SeqIO
from Bio.Seq import Seq
from operon_analyzer import genes


def load_sequence(operon: genes.Operon) -> Optional[Seq]:
    """ Loads the DNA sequence for a given operon's contig. """
    with gzip.open(operon.contig_filename, 'rt') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            if record.id == operon.contig:
                return Seq(str(record.seq))
