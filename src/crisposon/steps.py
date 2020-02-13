from Bio.Blast.Applications import NcbiblastpCommandline, NcbipsiblastCommandline

from crisposon.utils import get_neighborhood_ranges
from crisposon.parsers import parse_blast_xml

import os

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
        blast_cline = NcbipsiblastCommandline(query=self.orfs, db=self.db, evalue=self.e_val, outfmt=5, max_target_seqs=1, out=blast_out, num_threads=4)
        blast_cline()

        return blast_out
    
    def execute(self, orfs):

        self.orfs = orfs
        blast_out = self.run_blast()
        self.hits = parse_blast_xml(blast_out, self.name)

class Blastp(BlastStep):

    def __init__(self, db, name, e_val, working_dir):

        BlastStep.__init__(self, db, name, e_val, working_dir)
    
    def run_blast(self):

        file_name = "{}_blast.xml".format(self.name)
        blast_out = os.path.join(self.working_dir, file_name)
        blast_cline = NcbiblastpCommandline(query=self.orfs, db=self.db, evalue=self.e_val, outfmt=5, max_target_seqs=1, out=blast_out, num_threads=4)
        blast_cline()

        return blast_out
    
    def execute(self, orfs):

        self.orfs = orfs
        blast_out = self.run_blast()
        self.hits = parse_blast_xml(blast_out, self.name)

class SeedBlastpsi(Blastpsi):

    def __init__(self, db, name, e_val, working_dir, span):

        Blastpsi.__init__(self, db, name, e_val, working_dir)
        self.span = span
        self.is_seed = True

    def execute(self, orfs):
        
        self.orfs = orfs
        blast_out = self.run_blast()
        self.hits = parse_blast_xml(blast_out, self.name)
        self.neighborhood_ranges = get_neighborhood_ranges(self.hits, self.span)

class SeedBlastp(Blastp):

    def __init__(self, db, name, e_val, working_dir, span):

        Blastp.__init__(self, db, name, e_val, working_dir)
        self.span = span
        self.is_seed = True
    
    def execute(self, orfs):
        
        self.orfs = orfs
        blast_out = self.run_blast()
        self.hits = parse_blast_xml(blast_out, self.name)
        self.neighborhood_ranges = get_neighborhood_ranges(self.hits, self.span)

class FilterBlastpsi(Blastpsi):

    def __init__(self, db, name, e_val, working_dir):

        Blastp.__init__(self, db, name, e_val, working_dir)
        self.is_filter = True

    def execute(self, orfs):
        
        self.orfs = orfs
        blast_out = self.run_blast()
        self.hits = parse_blast_xml(blast_out, self.name)

class FilterBlastp(Blastp):

    def __init__(self, db, name, e_val, working_dir):

        Blastp.__init__(self, db, name, e_val, working_dir)
        self.is_filter = True

    def execute(self, orfs):
        
        self.orfs = orfs
        blast_out = self.run_blast()
        self.hits = parse_blast_xml(blast_out, self.name)
