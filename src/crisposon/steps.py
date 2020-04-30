from Bio.Blast.Applications import NcbiblastpCommandline, NcbipsiblastCommandline

from crisposon.utils import get_neighborhood_ranges
from crisposon.parsers import parse_blast_xml, parse_pilercr_summary, parse_mmseqs

import os, subprocess, tempfile, multiprocessing

class Blastp():
    """Wrapper for NCBI Blastp command line util.

    Really is just a wrapper around biopython's
    better wrapper.
    """

    def __init__(self, db, e_val, step_id):

        self.tmp_dir = tempfile.TemporaryDirectory()
        self.db = db
        self.e_val = e_val
        self.step_id = step_id

    def run(self, orfs):

        cores = multiprocessing.cpu_count()
        blast_out = os.path.join(self.tmp_dir.name, "blast_out.xml")
        blast_cline = NcbiblastpCommandline(query=orfs, db=self.db, evalue=self.e_val, 
                                                    outfmt=5, max_target_seqs=1, out=blast_out,
                                                    num_threads=cores)

        blast_cline()
        hits = parse_blast_xml(blast_out, self.step_id)
        return hits
        
class Blastpsi():
    """Wrapper for NCBI psiBlast command line util.

    Really is just a wrapper around biopython's
    better wrapper.
    """

    def __init__(self, db, e_val, step_id):

        self.tmp_dir = tempfile.TemporaryDirectory()
        self.db = db
        self.e_val = e_val
        self.step_id = step_id

    def run(self, orfs):

        cores = multiprocessing.cpu_count()
        blast_out = os.path.join(self.tmp_dir.name, "blast_out.xml")
        blast_cline = NcbipsiblastCommandline(query=orfs, db=self.db, evalue=self.e_val, 
                                                    outfmt=5, max_target_seqs=1, out=blast_out,
                                                    num_threads=cores)
        
        blast_cline()
        hits = parse_blast_xml(blast_out, self.step_id)
        return hits

class MMseqs():
    """Wrapper for mmseqs command line search util."""

    def __init__(self, db, e_val, name, sensitivity, extra_args=None):
        """Initialize a mmseqs command line run"""

        self.tmp_dir = tempfile.TemporaryDirectory()
        self.target_db = db
        self.name = name
        self.e_val = e_val
        self.sensitivity = sensitivity
        self.extra_args = extra_args
    
    def _make_query_db(self, orfs):
        """Covert query orfs to mmseqs custom database format."""

        query_db = "{}/querydb".format(self.tmp_dir.name)
        createdb_cmd = ["mmseqs", "createdb", orfs, query_db]
        subprocess.run(createdb_cmd, check=True)
        return query_db
    
    def _extract_best_hits(self, query_db, result_db):
        """Extract best hit info from raw mmseqs output.
        
        For simplicity, "best" means the hit with the 
        lowest e-value score.
        """

        results_tsv = "{}/results.tsv".format(self.tmp_dir.name)
        tsv_fields = "query,target,evalue,qseq,qheader,theader,qcov,tset"

        # convert raw output from custom mmseqs format to tsv format
        convertalis_cmd = ["mmseqs", "convertalis", query_db, self.target_db, 
                            result_db, results_tsv, "--format-output", tsv_fields]
        subprocess.run(convertalis_cmd, check=True)
        
        hits = parse_mmseqs(results_tsv, self.name, tsv_fields.split(","))
        return hits
    
    def run(self, orfs):

        # setup necessary files and directories for mmseqs run
        query_db = self._make_query_db(orfs)
        result_db = "{}/resultdb".format(self.tmp_dir.name)
        mmseqs_work = os.path.join(self.tmp_dir.name, "tmp")
        os.mkdir(mmseqs_work)

        cmd = ["mmseqs", "search", query_db, self.target_db, result_db, mmseqs_work, "-s",
                self.sensitivity, "-e", self.e_val, "-v", "0"]
        
        if self.extra_args is not None:
            cmd = cmd + self.extra_args
        subprocess.run(cmd, check=True)
        
        hits = self._extract_best_hits(query_db, result_db)
        return hits
    
class Diamond():
    """Wrapper for diamond command line search util."""

    def __init__(self, db, e_val, name, sensitivity, extra_args=None):
        """Initialize a mmseqs command line run"""

        self.tmp_dir = tempfile.TemporaryDirectory()
        self.db = db
        self.name = name
        self.e_val = e_val
        self.sensitivity = sensitivity
        self.extra_args = extra_args

    def run(self, orfs):

        result = os.path.join(self.tmp_dir.name, "result.xml")
        cmd = ["diamond", "blastp", "-q", orfs, "-d", self.db, "-o", result,
                self.sensitivity, "-e", self.e_val, "-k", "1", "-f", "5", "--quiet"]
        
        if self.extra_args is not None:
            cmd = cmd + self.extra_args
        subprocess.run(cmd, check=True)
        
        hits = parse_blast_xml(result, self.name)
        return hits       

class SearchStep:

    def __init__ (self, working_dir, search_tool):

        self.working_dir = working_dir
        self.search_tool = search_tool

    def execute(self, orfs):

        self.orfs = orfs
        self.hits = self.search_tool.run(orfs)

class SeedStep(SearchStep):

    def __init__ (self, working_dir, search_tool):

        SearchStep.__init__(self, working_dir, search_tool)
    
    def execute(self, orfs, span):

        super().execute(orfs)

        if len(self.hits) != 0:
            self.neighborhood_ranges = get_neighborhood_ranges(self.hits, span)
        else:
            self.neighborhood_ranges = None

class FilterStep(SearchStep):
    
    def __init__ (self, working_dir, search_tool, min_prot_count):

        SearchStep.__init__(self, working_dir, search_tool)
        self.min_prot_count = min_prot_count

class CrisprStep:

    def __init__(self, working_dir):

        self.working_dir = working_dir
        self.genome = None

    def run_pilercr(self):

        pilercr_out = os.path.join(self.working_dir, "pilercr_results")
        cmd = ["pilercr", "-in", self.genome, "-out", pilercr_out, "-minarray", "2", "-quiet"]
        subprocess.run(cmd, check=True)
        return pilercr_out
    
    def execute(self, genome):

        self.genome = genome
        pilercr_out = self.run_pilercr()
        self.hits = parse_pilercr_summary(pilercr_out)
