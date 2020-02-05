from crisposon.orffinder import orffinder, neighborhood_orffinder
from crisposon.utils import concatenate
from crisposon.dbbuilders import build_blastp_db
from crisposon.steps import SeedBlastp, SeedBlastpsi, FilterBlastp, FilterBlastpsi, Blastp, Blastpsi

import tempfile, os

BLASTP_KEYWORDS = ['Blastp', 'BlastP', 'BLASTP', 'PROT']
PSIBLAST_KEYWORDS = ['psiBlast', 'psiblast', 'PSIBLAST', 'PSI']

class Pipeline:

    def __init__(self, genome, id, outdir=None, min_prot_len=30, span=20000):

        self.genome = genome
        self.id = id
        self.outdir = outdir
        self.min_prot_len = min_prot_len
        self.span = span
        
        self._steps = []
        self._working_dir = tempfile.TemporaryDirectory()
        self._neighborhood_orfs = {}
        self._results = {}
        self._get_all_orfs()

    def __del__(self):

        self._working_dir.cleanup()
    
    def _results_init(self, neighborhood_ranges):

        for r in neighborhood_ranges:
            key = "{}_{}".format(r[0], r[1])
            self._results[key] = {"n_start": int(r[0]), "n_stop": int(r[1]), "marked": False, "hits": {}}
    
    def _results_update(self, hits):

        for hit in hits.keys():
            for neighborhood in self._results.keys():
                h_start = min(int(hits[hit]["q_start"]), int(hits[hit]["q_stop"]))
                h_stop = max(int(hits[hit]["q_start"]), int(hits[hit]["q_stop"]))
                if h_start >= self._results[neighborhood]["n_start"] and h_stop <= self._results[neighborhood]["n_stop"]:
                    self._results[neighborhood]["hits"][hit] = hits[hit]
        
    def _results_filter(self, hits):

        for hit in hits.keys():
            for neighborhood in self._results.keys():
                h_start = min(int(hits[hit]["q_start"]), int(hits[hit]["q_stop"]))
                h_stop = max(int(hits[hit]["q_start"]), int(hits[hit]["q_stop"]))
                if h_start >= self._results[neighborhood]["n_start"] and h_stop <= self._results[neighborhood]["n_stop"]:
                    self._results[neighborhood]["hits"][hit] = hits[hit]
                    self._results[neighborhood]["marked"] = True
        
        remove = [neighborhood for neighborhood in self._results if not self._results[neighborhood]["marked"]]
        for neighborhood in remove:
            del self._results[neighborhood]
            del self._neighborhood_orfs[neighborhood]
        
        for neighborhood in self._results.keys():
            self._results[neighborhood]["marked"] = False

    def _get_all_orfs(self):

        self._all_orfs = os.path.join(self._working_dir.name, "all_orfs.fasta")
        orffinder(sequence=self.genome, output=self._all_orfs, min_prot_len=self.min_prot_len, description=self.id)

    def _get_orfs_in_neighborhood(self, ranges):

        neighborhood_orffinder(sequence=self.genome, ranges=ranges, outdir=self._working_dir.name, min_prot_len=self.min_prot_len, description=self.id)

        for r in ranges:
            key = "{}_{}".format(r[0], r[1])
            path = os.path.join(self._working_dir.name, "orf_{}_{}.fasta".format(r[0], r[1]))
            self._neighborhood_orfs[key] = path

    def add_blast_seed_step(self, db, name, e_val, blast_type):

        if blast_type in BLASTP_KEYWORDS:
            self._steps.append(SeedBlastp(db, name, e_val, self._working_dir.name, self.span))
        elif blast_type in PSIBLAST_KEYWORDS:
            self._steps.append(SeedBlastpsi(db, name, e_val, self._working_dir.name, self.span))
        else:
            raise ValueError("blast type option '{}' not available for seed step".format(blast_type))

    def add_blast_filter_step(self, db, name, e_val, blast_type):

        if blast_type in BLASTP_KEYWORDS:
            self._steps.append(FilterBlastp(db, name, e_val, self._working_dir.name))
        elif blast_type in PSIBLAST_KEYWORDS:
            self._steps.append(FilterBlastpsi(db, name, e_val, self._working_dir.name))
        else:
            raise ValueError("blast type option '{}' not available for filter step".format(blast_type))
    
    def add_blast_step(self, db, name, e_val, blast_type):

        if blast_type in BLASTP_KEYWORDS:
            self._steps.append(Blastp(db, name, e_val, self._working_dir.name))
        elif blast_type in PSIBLAST_KEYWORDS:
            self._steps.append(Blastpsi(db, name, e_val, self._working_dir.name))
        else:
            raise ValueError("blast type option '{}' not available for filter step".format(blast_type))

    def run(self):

        neighborhood_orfs = None
        for step in self._steps:
            
            if step.is_seed:
                step.execute(self._all_orfs)
                self._get_orfs_in_neighborhood(step.neighborhood_ranges)
                self._results_init(step.neighborhood_ranges)
                self._results_update(step.hits)
                neighborhood_orfs = concatenate(self._working_dir.name, self._neighborhood_orfs.values())

            elif step.is_filter:
                step.execute(neighborhood_orfs)
                self._results_filter(step.hits)
                neighborhood_orfs = concatenate(self._working_dir.name, self._neighborhood_orfs.values())
            
            else:
                step.execute(neighborhood_orfs)
                self._results_update(step.hits)

        results = self._results
        return results


if __name__ == "__main__":

    import json

    genome = "/home/alexis/Projects/CRISPR-Transposons/data/genomes/v_crass_J520_whole.fasta"
    out = "/home/alexis/Projects/CRISPR-Transposons/data/"
    seed_db = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/tnsAB/blast_db"
    filter_db = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/cas_uniref/blast_db"
    final_db = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/tns_dc/blast_db"

    """in_dir = "/home/alexis/Projects/CRISPR-Transposons/data/protein_references/tns_cd/"
    db_dir = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/tns_dc"
    build_blastp_db(input=in_dir, db_dir=db_dir)"""
    
    p = Pipeline(genome, "v_crass", out, span=10000)
    p.add_blast_seed_step(seed_db, "tnsAB", 0.001, "PSI")
    p.add_blast_filter_step(filter_db, "cas", 0.001, "PROT")
    p.add_blast_step(final_db, "tnsCD", 0.001, "PROT")
    results = p.run()

    print(json.dumps(results, indent=4))
