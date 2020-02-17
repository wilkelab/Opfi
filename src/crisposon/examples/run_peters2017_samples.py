from crisposon.pipeline import Pipeline
import json, os

if __name__ == "__main__":
    
    # paths should probably be in a config file
    seed_db = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/cas6_tniQ/blast_db"
    filter_db1 = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/tnsAB/blast_db"
    filter_db2 = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/cas578/blast_db"
    tnsCD_db = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/tns_dc/blast_db"
    
    # genome sequences from the 14 organisms in fig3 peters et al 2017
    org_dir = "/home/alexis/Projects/CRISPR-Transposons/data/genomes/peters_2017_fig3"
    
    results_final = {}
    for genome in os.listdir(org_dir):

        path_to_infile = os.path.join(org_dir, genome)
        name = genome.strip(".fasta")
        results_final[name] = {}

        p = Pipeline(path_to_infile, name, min_prot_len=30, span=10000)
        p.add_seed_step(seed_db, "cas6_tniQ", 0.001, "PROT")
        p.add_filter_step(filter_db1, "tnsAB", 0.001, "PROT")
        p.add_filter_step(filter_db2, "cas578", 0.001, "PSI")
        p.add_blast_step(tnsCD_db, "tnsCD", 0.001, "PROT")
        p.add_crispr_step()
        results = p.run()

        for neighborhood in results.keys():
            results_final[name][neighborhood] = results[neighborhood]

    result = "/home/alexis/Projects/CRISPR-Transposons/data/misc/peters2017_replicated_results/cas6tniQ_seed"
    with open(result, 'w') as outfile:
        json.dump(results_final, outfile)