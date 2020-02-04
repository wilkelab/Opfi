from crisposon.orffinder import orffinder, neighborhood_orffinder
from crisposon.blastutils import get_neighborhood_ranges, concatenate
from crisposon.xmlparser import parse_blast

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline, NcbipsiblastCommandline

import tempfile, os

BLASTP_KEYWORDS = ['Blastp', 'BlastP', 'BLASTP', 'PROT']
PSIBLAST_KEYWORDS = ['psiBlast', 'psiblast', 'PSIBLAST', 'PSI']

class BlastStep:

    def __init__ (self, db, name, e_val, working_dir):

        self.db = db
        self.name = name
        self.e_val = e_val
        self.working_dir = working_dir
        self.orfs = None
        self.is_seed = False
        self.is_filter = False

class Blastpsi(BlastStep):

    def __init__(self, db, name, e_val, working_dir):

        BlastStep.__init__(self, db, name, e_val, working_dir)
    
    def run_blast(self):

        file_name = "{}_blast.xml".format(self.name)
        blast_out = os.path.join(self.working_dir, file_name)
        blast_cline = NcbipsiblastCommandline(query=self.orfs, db=self.db, evalue=self.e_val, outfmt=5, max_target_seqs=1, out=blast_out)
        blast_cline()

        return blast_out

class Blastp(BlastStep):

    def __init__(self, db, name, e_val, working_dir):

        BlastStep.__init__(self, db, name, e_val, working_dir)
    
    def run_blast(self):

        file_name = "{}_blast.xml".format(self.name)
        blast_out = os.path.join(self.working_dir, file_name)
        blast_cline = NcbiblastpCommandline(query=self.orfs, db=self.db, evalue=self.e_val, outfmt=5, max_target_seqs=1, out=blast_out)
        blast_cline()

        return blast_out

class SeedBlastpsi(Blastpsi):

    def __init__(self, db, name, e_val, working_dir, span):

        Blastpsi.__init__(self, db, name, e_val, working_dir)
        self.span = span
        self.is_seed = True

    def execute(self, orfs):
        
        self.orfs = orfs
        blast_out = self.run_blast()
        self.hits = parse_blast(blast_out, self.name)
        self.neighborhood_ranges = get_neighborhood_ranges(self.hits, self.span)

class SeedBlastp(Blastp):

    def __init__(self, db, name, e_val, working_dir, span):

        Blastp.__init__(self, db, name, e_val, working_dir)
        self.span = span
        self.is_seed = True
    
    def execute(self, orfs):
        
        self.orfs = orfs
        blast_out = self.run_blast()
        self.hits = parse_blast(blast_out, self.name)
        self.neighborhood_ranges = get_neighborhood_ranges(self.hits, self.span)

class FilterBlastpsi(Blastpsi):

    def __init__(self, db, name, e_val, working_dir):

        Blastp.__init__(self, db, name, e_val, working_dir)
        self.is_filter = True

    def execute(self, orfs):
        
        self.orfs = orfs
        blast_out = self.run_blast()
        self.hits = parse_blast(blast_out, self.name)

class FilterBlastp(Blastp):

    def __init__(self, db, name, e_val, working_dir):

        Blastp.__init__(self, db, name, e_val, working_dir)
        self.is_filter = True

    def execute(self, orfs):
        
        self.orfs = orfs
        blast_out = self.run_blast()
        self.hits = parse_blast(blast_out, self.name)

class Pipeline:

    def __init__(self, genome, name, outdir, min_prot_len=30, span=20000):

        self.genome = genome
        self.name = name
        self.outdir = outdir
        self.min_prot_len = min_prot_len
        self.span = span
        self._working_dir = tempfile.TemporaryDirectory()

        self._neighborhood_orfs = {}
        self._steps = []
        self.results = {}
        self._get_all_orfs()

    def __del__(self):

        self._working_dir.cleanup()
    
    def _results_init(self, ranges):

        for rang in ranges:
            key = "{}_{}".format(rang[0], rang[1])
            self.results[key] = {"neighborhood_start": int(rang[0]), "neighborhood_stop": int(rang[1]), "keep": False, "hits": {}}
    
    def _results_update(self, hits):

        for hit in hits.keys():
            for neighborhood in self.results.keys():
                h_start = min(int(hits[hit]["query_start"]), int(hits[hit]["query_stop"]))
                h_stop = max(int(hits[hit]["query_start"]), int(hits[hit]["query_stop"]))
                if h_start >= self.results[neighborhood]["neighborhood_start"] and h_stop <= self.results[neighborhood]["neighborhood_stop"]:
                    self.results[neighborhood]["hits"][hit] = hits[hit]
        
    def _results_filter(self, hits):

        for hit in hits.keys():
            for neighborhood in self.results.keys():
                h_start = min(int(hits[hit]["query_start"]), int(hits[hit]["query_stop"]))
                h_stop = max(int(hits[hit]["query_start"]), int(hits[hit]["query_stop"]))
                if h_start >= self.results[neighborhood]["neighborhood_start"] and h_stop <= self.results[neighborhood]["neighborhood_stop"]:
                    self.results[neighborhood]["hits"][hit] = hits[hit]
                    self.results[neighborhood]["keep"] = True
        
        delete = [neighborhood for neighborhood in self.results if not self.results[neighborhood]["keep"]]
        for neighborhood in delete:
            del self.results[neighborhood]
            del self._neighborhood_orfs[neighborhood]
        
        for neighborhood in self.results.keys():
            self.results[neighborhood]["keep"] = False

    def _get_all_orfs(self):

        self._all_orfs = os.path.join(self._working_dir.name, "all_orfs.fasta")
        orffinder(sequence=self.genome, output=self._all_orfs, min_prot_len=self.min_prot_len, description=self.name)

    def _get_neighborhood_orfs(self, ranges):

        neighborhood_orffinder(sequence=self.genome, ranges=ranges, outdir=self._working_dir.name, min_prot_len=self.min_prot_len, description=self.name)

        for rang in ranges:
            key = "{}_{}".format(rang[0], rang[1])
            path = os.path.join(self._working_dir.name, "orf_{}_{}.fasta".format(rang[0], rang[1]))
            self._neighborhood_orfs[key] = path

    def add_seed_step(self, db, name, e_val, blast_type):

        if blast_type in BLASTP_KEYWORDS:
            self._steps.append(SeedBlastp(db, name, e_val, self._working_dir.name, self.span))
        elif blast_type in PSIBLAST_KEYWORDS:
            self._steps.append(SeedBlastpsi(db, name, e_val, self._working_dir.name, self.span))
        else:
            raise ValueError("blast type option '{}' not available for seed step".format(blast_type))

    def add_filter_step(self, db, name, e_val, blast_type):

        if blast_type in BLASTP_KEYWORDS:
            self._steps.append(FilterBlastp(db, name, e_val, self._working_dir.name))
        elif blast_type in PSIBLAST_KEYWORDS:
            self._steps.append(FilterBlastpsi(db, name, e_val, self._working_dir.name))
        else:
            raise ValueError("blast type option '{}' not available for filter step".format(blast_type))

    def run(self):

        for step in self._steps:
            
            if len(self._neighborhood_orfs) > 0:
                neighborhood_orfs = concatenate(self._working_dir.name, self._neighborhood_orfs.values())
            else:
                neighborhood_orfs = None
            
            if step.is_seed:
                step.execute(self._all_orfs)
                self._get_neighborhood_orfs(step.neighborhood_ranges)
                self._results_init(step.neighborhood_ranges)
                self._results_update(step.hits)

            elif step.is_filter:
                step.execute(neighborhood_orfs)
                self._results_filter(step.hits)

        results = self.results
        return results


if __name__ == "__main__":

    import json

    genome = "/home/alexis/Projects/CRISPR-Transposons/data/genomes/v_crass_J520_whole.fasta"
    out = "/home/alexis/Projects/CRISPR-Transposons/data/"
    seed_db = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/tnsAB/blast_db"
    filter_db = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/cas_uniref/blast_db"
    
    p = Pipeline(genome, "v_crass", out, span=5000)
    p.add_seed_step(seed_db, "tnsAB", 0.001, "PSI")
    p.add_filter_step(filter_db, "cas", 0.001, "PROT")
    results = p.run()

    print(json.dumps(results, indent=4))