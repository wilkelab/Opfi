from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def build_protein_fasta(orfs, out_dir, description):
    """
    Write orfs to a single fasta file.
    """
    records = []
    # sort orfs by position in parent
    orfs.sort(key=lambda x: min(x[0][0], x[0][1]))
    
    for orf in orfs:
        prot_id = "lcl|" + "|".join(str(x) for x in orf[0])
        prot_seq = Seq(orf[1])
        records.append(SeqRecord(prot_seq, id=prot_id, description=description))

    SeqIO.write(records, out_dir, "fasta")


def aa_index_conversion(aa_index, frame, strand, parent_len, start=True):
    """
    Convert the index of an amino acid to it's position in the parent 
    nucleotide sequence.
    """
    nuc_index = 0
    
    if strand == 1:
        if start:
            nuc_index = aa_index * 3 + 1 + frame
        else: 
            nuc_index = aa_index * 3 + 3 + frame
    else:
        if start:
            nuc_index = parent_len - aa_index * 3 - frame
        else:
            nuc_index = parent_len - aa_index * 3 - 2 - frame
    
    return nuc_index


def get_orfs_in_range(all_orfs, rang):
    """
    Given a list of orfs, return only those contained within a range.
    """
    orfs_in_range = []
    for orf in all_orfs:
        lower = min(orf[0][0], orf[0][1])
        upper = max(orf[0][0], orf[0][1])

        if lower >= rang[0] and upper <= rang[1]:
            orfs_in_range.append(orf)

    return orfs_in_range


def get_orfs_in_frame(seq, frame, strand, min_prot_len, all_orfs):
    """
    Get all orfs in a single reading frame.
    """
    parent_len = len(seq)
    # Add extra Ns to the sequence to ensure we have a complete frame
    buffer_len = 3 - ((parent_len - frame) % 3)
    sequence_buffer = "NNN"[:buffer_len]

    trans = str((seq[frame:] + sequence_buffer).translate())
    trans_len = len(trans)
    aa_start = 0

    while aa_start < trans_len:
        # find the next stop site
        aa_stop = trans.find("*", aa_start)
        if aa_stop == -1: # no stop sites remain
            aa_start = trans_len
        else:
            # find the methionine that makes this ORF as long as possible
            aa_start = trans.find("M", aa_start, aa_stop)
            prot_len = aa_stop - aa_start

            # if a methionine was found and the ORF is longer than the cutoff, 
            # store it's data in 'all_orfs'
            if aa_start != -1 and prot_len >= min_prot_len:
                nuc_start = aa_index_conversion(aa_start, frame, strand, parent_len, start=True)
                nuc_end = aa_index_conversion(aa_stop, frame, strand, parent_len, start=False)
                
                prot_seq = trans[aa_start:aa_stop]
                orf_data = [nuc_start, nuc_end, frame + 1, strand]
                orf = (orf_data, prot_seq)
                all_orfs.append(orf)
            aa_start = aa_stop + 1


def get_all_orfs(seq, min_prot_len):
    """
    Get orfs in all six reading frames.
    """
    all_orfs = []

    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            get_orfs_in_frame(nuc, frame, strand, min_prot_len, all_orfs)
    return all_orfs


def reader(contig):
    try:
        record = SeqIO.read(contig, "fasta")
    except:
        print("Cannot read input file. Input should be in fasta format and contain only one record.")
    return record.seq


def orffinder(sequence, output, min_prot_len=30, description="putative protein"):
    """
    Find all open reading frames in a nucleotide sequence.

    ORFs are translated to amino acid sequences and written to a protein
    fasta. The unique identifier for each translated ORF contains information
    about it's start and end positions in the parent sequence, the reading 
    frame (1, 2, or 3), and the strand (1 for original sequence, -1 for
    the reverse compliment). ID format: lcl|<start>|<end>|<frame>|<strand> 

    Args:
        sequence (str): Path to the input sequence file (in fasta format).
        output (str): Path to the output file.
        min_prot_len (int, optional): Minimal ORF length (amino acids). 
            Defaults to 30 residues.
        description (str, optional): Will be applied to all output records. 
            Defaults to "putative protein".
    """
    parent_seq = reader(sequence)
    all_orfs = get_all_orfs(parent_seq, min_prot_len)
    build_protein_fasta(all_orfs, output, description)
    if len(all_orfs) == 0:
        return None
    else:
        return output


def neighborhood_orffinder(sequence, ranges, outdir, min_prot_len=30, description="putative protein"):
    """
    Find open reading frames in a nucleotide sequence that occur within a certain range/ranges.

    Args:
        sequence (str): Path to the input sequence file (in fasta format).
        ranges (list): List of tuples containing the start and stop position
            for each range.
        outdir (str): Path to the output directory.
        min_prot_len (int, optional): Minimal ORF length (amino acids). 
            Defaults to 30 residues.
        description (str, optional): Will be applied to all output records. 
            Defaults to "putative protein".
    """
    parent_seq = reader(sequence)
    all_orfs = get_all_orfs(parent_seq, min_prot_len)

    for rang in ranges:
        orfs_in_range = get_orfs_in_range(all_orfs, rang)
        outfil = "orf_{}_{}.fasta".format(str(rang[0]), str(rang[1]))
        out = os.path.join(outdir, outfil)
        build_protein_fasta(orfs_in_range, out, description)
