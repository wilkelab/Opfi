from crisposon.pipeline import Pipeline
from Bio import SeqIO
import yaml
import json, os, sys, tempfile

if __name__ == "__main__":
    """Example usage of the crisposon package

    Command line args:
        1. Path to input fasta file
        2. Path to steps config file
    
    Constructs a Pipeline object to search for CRISPR-Transposon 
    systems in the genome/config of interest.
    
    """
    
    genomic_data = sys.argv[1]
    conf_file = sys.argv[2]

    # load contents of yaml conf file
    stream = open(conf_file, 'r')
    conf = yaml.load(stream)
    
    work = tempfile.TemporaryDirectory()
    
    job_id = os.path.basename(genomic_data)
    job_out = os.path.join(conf["out-dir"], job_id)
    if not os.path.exists(job_out):
        os.mkdir(job_out)

    records = SeqIO.parse(genomic_data, "fasta")
    for record in records:

        # each record needs to be written to its own fasta file
        # will be automatically deleted when the program finishes
        contig = os.path.join(work.name, "{}.fasta".format(record.id))
        SeqIO.write(record, contig, "fasta")

        p = Pipeline(genome=contig, id=record.id, min_prot_len=conf["min-prot-len"], span=conf["span"])
        for step in conf["steps"]:
            
            if step["type"] == "seed":
                p.add_seed_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"])
            
            elif step["type"] == "filter":
                p.add_filter_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"])
            
            elif step["type"] == "blast":
                p.add_blast_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"])
            
            else:
                p.add_crispr_step()
        
        results = p.run()
        
        if len(results) != 0:
            result_out = os.path.join(job_out, "{}.json".format(record.id))
            with open(result_out, 'w') as outfile:
                json.dump(results, outfile)
