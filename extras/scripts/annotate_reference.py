from Bio import SeqIO
import os, sys

def annotate_reference(prot_ref_file, label):
    """
    Prepends each record description in the protein reference file with 
    a label (the gene name).

    Makes downstream processing/visualization easier.
    """
    records = list(SeqIO.parse(ref_fasta, "fasta"))
        
    for record in records:
        des = record.description.split()
        prot_id = des.pop(0)
        des_with_label = "{} {} {}".format(prot_id, label, " ".join(des))
        record.description = des_with_label

    SeqIO.write(records, ref_fasta, "fasta")

if __name__ == "__main__":
    ref_fasta = sys.argv[1]
    label = sys.argv[2]
    annotate_reference(ref_fasta, label)