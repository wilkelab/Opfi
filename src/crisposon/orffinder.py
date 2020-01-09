from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def build_protein_fasta(orfs, out_dir, description):
    records = []
    orfs.sort(key=lambda x: len(x[1]), reverse=True) # sort by seq len in desc order
    for orf in orfs:
        prot_id = "lcl|" + "|".join(str(x) for x in orf[0])
        prot_seq = Seq(orf[1])
        records.append(SeqRecord(prot_seq, id=prot_id, description=description))

    SeqIO.write(records, out_dir, "fasta")

# Convert the index of an amino acid to it's position in the parent 
# nucleotide sequence.
def aa_index_conversion(aa_index, frame, strand, parent_len, start=True):
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

def get_orfs_in_frame(seq, frame, strand, min_prot_len, all_orfs):
    parent_len = len(seq)
    trans = str(seq[frame:].translate())
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
    all_orfs = []

    # Look for orfs in all six reading frames
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

# Actual public function to call
def orffinder(sequence, output, min_prot_len=30, description="putative protein"):
    """Find all open reading frames in a nucleotide sequence.

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