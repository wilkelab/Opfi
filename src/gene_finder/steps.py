import multiprocessing
import os
import subprocess
import tempfile

from Bio.Blast.Applications import (NcbiblastpCommandline,
                                    NcbipsiblastCommandline)

from gene_finder.option_handling import (build_blastn_command,
                                         build_blastp_command,
                                         build_psiblast_command)
from gene_finder.parsers import (parse_blastn_output, parse_pilercr_summary,
                                 parse_search_output)
from gene_finder.utils import get_neighborhood_ranges


class Blastp():
    """
    Wrapper for the NCBI BLASTP command line utility.
    """
    BLASTOUT_FIELDS = "qseqid sseqid stitle evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovhsp qseq"

    def __init__(self, db, e_val, step_id, parse_descriptions, blastp_path, kwargs):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.db = db
        self.e_val = e_val
        self.step_id = step_id
        self.parse_descriptions = parse_descriptions
        self.blastp_path = blastp_path
        self.kwargs = kwargs
    
    def construct_cmd(self, query, out):
        """
        Format the BLASTP command into a string that `subprocess.run()` can use.
        """
        return build_blastp_command(query, self.db, self.e_val, self.kwargs, self.BLASTOUT_FIELDS, out, self.blastp_path)

    def run(self, orfs):
        """
        Execute the BLASTP search.
        """
        blast_out = os.path.join(self.tmp_dir.name, "blast_out.tsv")
        cmd = self.construct_cmd(orfs, blast_out)
        subprocess.run(cmd, check=True)
        hits = parse_search_output(blast_out, self.step_id, "blast", self.parse_descriptions)
        return hits
        

class Blastpsi():
    """
    Wrapper for NCBI psiBLAST command line utility.
    """
    BLASTOUT_FIELDS = "qseqid sseqid stitle evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovhsp qseq"

    def __init__(self, db, e_val, step_id, parse_descriptions, psiblast_path, kwargs):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.db = db
        self.e_val = e_val
        self.step_id = step_id
        self.parse_descriptions = parse_descriptions
        self.psiblast_path = psiblast_path
        self.kwargs = kwargs

    def construct_cmd(self, query, out):
        """
        Format the psiBLAST command into a string that `subprocess.run()` can use.
        """
        return build_psiblast_command(query, self.db, self.e_val, self.kwargs, self.BLASTOUT_FIELDS, out, self.psiblast_path)

    def run(self, orfs):
        """
        Execute the psiBLAST search.
        """
        blast_out = os.path.join(self.tmp_dir.name, "blast_out.tsv")
        cmd = self.construct_cmd(orfs, blast_out)
        subprocess.run(cmd, check=True)
        hits = parse_search_output(blast_out, self.step_id, "blast", self.parse_descriptions)
        return hits


class MMseqs():
    """
    Wrapper for mmseqs command line search utility.
    """
    MMSEQSOUT_FIELDS = "qheader target theader evalue bits raw alnlen pident nident mismatch gapopen qcov qseq"

    def __init__(self, db, e_val, step_id, sensitivity, parse_descriptions):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.target_db = db
        self.step_id = step_id
        self.e_val = e_val
        self.sensitivity = sensitivity
        self.parse_descriptions = parse_descriptions
    
    def _make_query_db(self, orfs):
        """
        Convert query orfs to mmseqs custom database format.
        """
        query_db = "{}/querydb".format(self.tmp_dir.name)
        createdb_cmd = ["mmseqs", "createdb", orfs, query_db]
        subprocess.run(createdb_cmd, check=True)
        return query_db
    
    def _extract_best_hits(self, query_db, result_db):
        """
        Extract best hit info from raw mmseqs output.
        For simplicity, "best" means the hit with the 
        lowest e-value score.
        """
        results_tsv = "{}/results.tsv".format(self.tmp_dir.name)
        tsv_fields = ",".join(self.MMSEQSOUT_FIELDS.split())

        # convert raw output from custom mmseqs format to tsv format
        convertalis_cmd = ["mmseqs", "convertalis", query_db, self.target_db, 
                            result_db, results_tsv, "--format-output", tsv_fields]
        subprocess.run(convertalis_cmd, check=True)
        hits = parse_search_output(results_tsv, self.step_id, "mmseqs", self.parse_descriptions)
        return hits
    
    def run(self, orfs):
        """
        Execute the mmseqs search and parse output.
        """
        # setup necessary files and directories for mmseqs run
        query_db = self._make_query_db(orfs)
        result_db = "{}/resultdb".format(self.tmp_dir.name)
        mmseqs_work = os.path.join(self.tmp_dir.name, "tmp")
        os.mkdir(mmseqs_work)

        cmd = ["mmseqs", "search", query_db, self.target_db, result_db, mmseqs_work, "-e", self.e_val, "-v", "0"]
        if self.sensitivity is not None:
            cmd.append("-s")
            cmd.append(str(self.sensitivity))
        subprocess.run(cmd, check=True)
        # parse output
        hits = self._extract_best_hits(query_db, result_db)
        return hits
    

class Diamond():
    """
    Wrapper for diamond command line search utility.
    """
    DIAMONDOUT_FIELDS = "qseqid sseqid stitle evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovhsp qseq"

    def __init__(self, db, e_val, step_id, sensitivity, parse_descriptions):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.db = db
        self.step_id = step_id
        self.e_val = e_val
        self.sensitivity = sensitivity
        self.parse_descriptions = parse_descriptions

    def run(self, orfs):
        """
        Execute the diamond blastp search and parse output.
        """
        result = os.path.join(self.tmp_dir.name, "result.tsv")
        cmd = ['diamond', 'blastp', '-q', orfs, '-d', self.db, '-o', result,
                '-e', self.e_val, '--threads', '1', '--quiet', '-f', '6']
        cmd = cmd + self.DIAMONDOUT_FIELDS.split()
        if self.sensitivity is not None:
            cmd.append(str(self.sensitivity))
        subprocess.run(cmd, check=True)
        # parse output
        hits = parse_search_output(result, self.step_id, "diamond", self.parse_descriptions)
        return hits       


class Pilercr():
    """
    Wrapper for PILER-CR CRISPR array finding software.
    """
    def __init__(self, step_id):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.step_id = step_id
    
    def run(self, genome):
        """
        Execute the PILER-CR array search. 
        """
        pilercr_out = os.path.join(self.tmp_dir.name, "pilercr_results")
        cmd = ["pilercr", "-in", genome, "-out", pilercr_out, "-minarray", "2", "-quiet"]
        subprocess.run(cmd, check=True)
        hits = parse_pilercr_summary(pilercr_out)
        return hits
        

class Blastn():
    """
    Wrapper for NCBI blastn command line utility.
    """
    BLASTOUT_FIELDS = "qseqid sseqid stitle evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovhsp qseq qstart qend sstrand"

    def __init__(self, db, step_id, e_val, parse_descriptions, blastn_path, kwargs):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.db = db
        self.step_id = step_id
        self.e_val = e_val
        self.parse_descriptions = parse_descriptions
        self.blastn_path = blastn_path
        self.kwargs = kwargs
    
    def construct_cmd(self, query, out):
        return build_blastn_command(query, self.db, self.e_val, self.kwargs, self.BLASTOUT_FIELDS, out, self.blastn_path)

    def run(self, genome):
        blast_out = os.path.join(self.tmp_dir.name, "blast_out.tsv")
        cmd = self.construct_cmd(genome, blast_out)
        subprocess.run(cmd, check=True)
        hits = parse_blastn_output(blast_out, self.step_id, self.parse_descriptions)
        return hits      


class SearchStep:
    """
    Represents a generic step in the pipeline.
    """
    def __init__ (self, search_tool):
        self.search_tool = search_tool

    def execute(self, orfs):
        self.orfs = orfs
        self.hits = self.search_tool.run(orfs)


class SeedStep(SearchStep):
    """
    Represents a bait/seed step in the pipeline.
    """
    def __init__ (self, search_tool):
        SearchStep.__init__(self, search_tool)
    
    def execute(self, orfs, span, contig_len):
        super().execute(orfs)
        if len(self.hits) != 0:
            self.neighborhood_ranges = get_neighborhood_ranges(self.hits, contig_len, span)
        else:
            self.neighborhood_ranges = []


class SeedWithCoordinatesStep():
    """
    Represents a bait/seed step in the pipeline, where coordinates are used to 
    define candidates.
    """
    def __init__(self, start, end, contig_id):
        # these shouldn't change after initialization
        self.start = start
        self.end = end
        self.contig_id = contig_id

        # Will be updated for each contig if no 
        # coordinates were given when the step was added
        self.coordinates = (start, end)
    
    def update_start_coord(self, new_start):
        self.coordinates = (new_start, self.coordinates[1])
    
    def update_end_coord(self, new_end):
        self.coordinates = (self.coordinates[0], new_end)
    
    @property
    def neighborhood_ranges(self):
        return [self.coordinates]


class FilterStep(SearchStep):
    """
    Represents a filtering step in the pipeline. 
    """
    def __init__ (self, search_tool, min_prot_count):
        SearchStep.__init__(self, search_tool)
        self.min_prot_count = min_prot_count


class CrisprStep(SearchStep):
    """
    Represents a CRISPR array search step in the pipeline.
    """
    def __init__(self, search_tool):
        SearchStep.__init__(self, search_tool)


class BlastnStep(SearchStep):
    """
    Represents a nucleotide BLAST search step in the pipeline.
    """
    def __init__(self, search_tool):
        SearchStep.__init__(self, search_tool)
