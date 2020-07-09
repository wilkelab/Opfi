from Bio.Blast.Applications import (NcbiblastpCommandline, 
                                    NcbipsiblastCommandline)

from gene_finder.utils import get_neighborhood_ranges
from gene_finder.parsers import parse_search_output, parse_pilercr_summary
from gene_finder.option_handling import (build_blastp_command,
                                        build_psiblast_command)

import os, subprocess, tempfile, multiprocessing

class Blastp():
    """Wrapper for NCBI Blastp command line util.

    Really is just a wrapper around biopython's
    better wrapper.
    """
    BLASTOUT_FIELDS = "qseqid sseqid stitle evalue \
        bitscore score length pident \
        nident mismatch positive gapopen \
        gaps ppos qcovhsp qseq"

    def __init__(self, db, e_val, step_id, kwargs):

        self.tmp_dir = tempfile.TemporaryDirectory()
        self.db = db
        self.e_val = e_val
        self.step_id = step_id
        self.kwargs = kwargs
    
    def construct_cmd(self, query, out):
        return build_blastp_command(query, self.db, self.e_val, self.kwargs, self.BLASTOUT_FIELDS, out)

    def run(self, orfs):

        blast_out = os.path.join(self.tmp_dir.name, "blast_out.tsv")
        cmd = self.construct_cmd(orfs, blast_out)
        subprocess.run(cmd, check=True)
        hits = parse_search_output(blast_out, self.step_id, "blast")
        return hits
        
class Blastpsi():
    """Wrapper for NCBI psiBlast command line util.

    Really is just a wrapper around biopython's
    better wrapper.
    """
    BLASTOUT_FIELDS = "qseqid sseqid stitle evalue \
        bitscore score length pident \
        nident mismatch positive gapopen \
        gaps ppos qcovhsp qseq"

    def __init__(self, db, e_val, step_id, kwargs):

        self.tmp_dir = tempfile.TemporaryDirectory()
        self.db = db
        self.e_val = e_val
        self.step_id = step_id
        self.kwargs = kwargs
    
    def construct_cmd(self, query, out):
        return build_psiblast_command(query, self.db, self.e_val, self.kwargs, self.BLASTOUT_FIELDS, out)

    def run(self, orfs):

        blast_out = os.path.join(self.tmp_dir.name, "blast_out.tsv")
        cmd = self.construct_cmd(orfs, blast_out)
        subprocess.run(cmd, check=True)
        hits = parse_search_output(blast_out, self.step_id, "blast")
        return hits

class MMseqs():
    """Wrapper for mmseqs command line search util."""

    MMSEQSOUT_FIELDS = "qheader target theader evalue \
        bits raw alnlen pident \
        nident mismatch gapopen qcov qseq"

    def __init__(self, db, e_val, step_id, sensitivity):
        """Initialize a mmseqs command line run"""

        self.tmp_dir = tempfile.TemporaryDirectory()
        self.target_db = db
        self.step_id = step_id
        self.e_val = e_val
        self.sensitivity = sensitivity
    
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
        tsv_fields = self.MMSEQSOUT_FIELDS.split().join(",")

        # convert raw output from custom mmseqs format to tsv format
        convertalis_cmd = ["mmseqs", "convertalis", query_db, self.target_db, 
                            result_db, results_tsv, "--format-output", tsv_fields]
        subprocess.run(convertalis_cmd, check=True)
        
        hits = parse_search_output(results_tsv, self.step_id, "mmseqs")
        return hits
    
    def run(self, orfs):
        """Execute the mmseqs search command and parse output."""

        # setup necessary files and directories for mmseqs run
        query_db = self._make_query_db(orfs)
        result_db = "{}/resultdb".format(self.tmp_dir.name)
        mmseqs_work = os.path.join(self.tmp_dir.name, "tmp")
        os.mkdir(mmseqs_work)

        cmd = ["mmseqs", "search", query_db, self.target_db, result_db, mmseqs_work, "-s",
                self.sensitivity, "-e", self.e_val, "-v", "0"]
        subprocess.run(cmd, check=True)
        
        # parse output
        hits = self._extract_best_hits(query_db, result_db)
        return hits
    
class Diamond():
    """Wrapper for diamond command line search util."""

    DIAMONDOUT_FIELDS = "qseqid sseqid stitle evalue \
        bitscore score length pident \
        nident mismatch positive gapopen \
        gaps ppos qcovhsp qseq"

    def __init__(self, db, e_val, step_id, sensitivity):
        """Initialize a diamond command line run"""

        self.tmp_dir = tempfile.TemporaryDirectory()
        self.db = db
        self.step_id = step_id
        self.e_val = e_val
        self.sensitivity = sensitivity

    def run(self, orfs):
        """Execute the diamond blastp command and parse output."""

        result = os.path.join(self.tmp_dir.name, "result.tsv")
        
        cmd = ["diamond", "blastp", "-q", orfs, "-d", self.db, "-o", result,
                    "-e", self.e_val, "-k", "1", "--quiet", "-f", 
                    "6", self.DIAMONDOUT_FIELDS]
        if self.sensitivity:
            cmd.append(self.sensitivity)
        subprocess.run(cmd, check=True, shell=True)
        
        # parse output
        hits = parse_search_output(result, self.step_id, "diamond")
        return hits       

class Pilercr():

    def __init__(self, step_id):

        self.tmp_dir = tempfile.TemporaryDirectory()
        self.step_id = step_id
    
    def run(self, genome):

        pilercr_out = os.path.join(self.tmp_dir.name, "pilercr_results")
        cmd = ["pilercr", "-in", genome, "-out", pilercr_out, "-minarray", "2", "-quiet"]
        subprocess.run(cmd, check=True)
        
        hits = parse_pilercr_summary(pilercr_out)
        return hits
        

class SearchStep:

    def __init__ (self, search_tool):

        self.search_tool = search_tool

    def execute(self, orfs):

        self.orfs = orfs
        self.hits = self.search_tool.run(orfs)

class SeedStep(SearchStep):

    def __init__ (self, search_tool):

        SearchStep.__init__(self, search_tool)
    
    def execute(self, orfs, span, contig_len):

        super().execute(orfs)

        if len(self.hits) != 0:
            self.neighborhood_ranges = get_neighborhood_ranges(self.hits, contig_len, span)
        else:
            self.neighborhood_ranges = None

class FilterStep(SearchStep):
    
    def __init__ (self, search_tool, min_prot_count):

        SearchStep.__init__(self, search_tool)
        self.min_prot_count = min_prot_count

class CrisprStep(SearchStep):

    def __init__(self, search_tool):

        SearchStep.__init__(self, search_tool)
