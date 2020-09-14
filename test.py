from gene_finder.pipeline import Pipeline


data = "tests/integration/integration_data/contigs/trna.fasta"
db = "tests/integration/integration_data/blast_databases/trna/trnas.fa"
p = Pipeline()
p.add_blastn_step(db, "test", 1e-3, False)
results = p.run(job_id="blast_test", data=data, output_directory="/tmp")
print(results)
