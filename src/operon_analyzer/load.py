import gzip
from Bio import SeqIO
from operon_analyzer import genes


def load_sequence(operon: genes.Operon):
    """ Loads an operon's nucleotide sequence from disk. """
    with gzip.open(operon.contig_filename, 'rt') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            if record.id == operon.contig:
                operon.set_sequence(record.seq)
                break
