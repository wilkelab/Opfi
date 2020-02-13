from crisposon.orffinder import orffinder, neighborhood_orffinder
from crisposon.utils import concatenate
from crisposon.build_blast_db import build_blastp_db
from crisposon.steps import SeedBlastp, SeedBlastpsi, FilterBlastp, FilterBlastpsi, Blastp, Blastpsi, CrisprStep

import tempfile, os

BLASTP_KEYWORDS = ['Blastp', 'BlastP', 'BLASTP', 'PROT']
PSIBLAST_KEYWORDS = ['psiBlast', 'psiblast', 'PSIBLAST', 'PSI']

class Pipeline:
    """
    Main class for running the CRISPR-transposon identification pipeline.
    
    Takes a single genome of interest as input, and performs a series of
     user-specified alignment steps to identify putative CRISPR-transposon elements.

    Args:
        genome (str): Path to genome/contig fasta file.
        id (str): Unique identifier for this genome/contig.
        min_prot_len (int, optional): Min residue size of ORFs
        in query genome. Default is 30.
        span (int, optional): Size of nt regions to keep around hits after the
        seed phase (see add_blast_seed_step method for details).

    Example:

        Create a pipeline to search for CRISPR-transposon systems
        in the Vibrio crassostreae genome.
        . 
        >>> p = Pipeline(genome="v_crass.fasta", id="v_crass", span=15000)

        Add alignment steps to the pipeline. First, we will do a blast
        for TnsA/B genes - "seeds" - in the query genome. Regions that fall outside
        of the span around each hit will be filtered out of subsequent searches.

        >>> p.add_blast_seed_step(db="blast_databases/tnsAB", name="tnsAB", e_val=0.001, blast_type="PROT")

        Now add a filter step for cas proteins. This will tell the program to
        to run a blast against a cas reference database, but only for the
        regions around the hits from the previous step. Regions ("neighborhoods")
        that do not contain hits for cas proteins will be filtered out of
        subsequent searches, as well as the results.

        >>> p.add_blast_filter_step(db="blast_databases/cas", name="cas", e_val=0.001, blast_type="PSI")

        Finally, add a blast step for tnsC and tnsD proteins using any remaining 
        neighborhoods as queries. Neighborhoods that do not contain hits will NOT
        be filtered out during this step.

        >>> p.add_blast_step(db="blast_databases/tnsCD", name="tnsCD", e_val=0.001, blast_type="PROT")

        Now, run the pipeline. Results are returned as a dictionary object containing
        the hits associated with each neighborhood. These are your putative CRISPR-transposon
        systems!

        >>> results = p.run()
    """

    def __init__(self, genome, id, min_prot_len=30, span=20000):
        """Initialize a Pipeline object with a genome/contig (path),
        unique id, and (optionally) a minimum ORF length and the 
        size of neighborhood regions.

        Initialize an empty list to keep track of steps (which will
        be added later) and an empty dictionary to track results.

        Also initialize an empty dictionary to keep track of the ORFs
        associated with each gene neighborhood. All blasts are
        executed using ORFs from the parent genome as queries, which
        helps to simplify filtering and results reporting.

        Set up a temporary working directory for intermediate files.
        """
        self.genome = genome
        self.id = id
        self.min_prot_len = min_prot_len
        self.span = span
        
        self._steps = []
        self._working_dir = tempfile.TemporaryDirectory()
        self._neighborhood_orfs = {}
        self._results = {}
        self._get_all_orfs()

    def __del__(self):
        """Delete the working directory and its contents when 
        this object is garbage collected.
        """
        self._working_dir.cleanup()
    
    def _results_init(self, neighborhood_ranges):
        """Create an entry for every gene neighborhood identified
        during the seed phase. All hits will be recorded here.
        """
        for r in neighborhood_ranges:
            key = "{}_{}".format(r[0], r[1])
            self._results[key] = {"n_start": int(r[0]), "n_stop": int(r[1]), "do_not_filter": False, "hits": {}}
    
    def _results_update(self, hits):
        """Add hits to the results tracker (grouped by neighborhood)."""

        for hit in hits.keys():
            for neighborhood in self._results.keys():

                # note that the begining and end of the hit is denoted by the begining and
                # end of the orf used as the query
                h_start = min(int(hits[hit]["Start"]), int(hits[hit]["Stop"]))
                h_stop = max(int(hits[hit]["Start"]), int(hits[hit]["Stop"]))
                
                # check whether this hit is contained within the neighborhood
                if h_start >= self._results[neighborhood]["n_start"] and h_stop <= self._results[neighborhood]["n_stop"]:
                    self._results[neighborhood]["hits"][hit] = hits[hit]
    
    def _results_update_crispr(self, hits):
        """Add hits for CRISPR arrays to the results tracker."""

        for hit in hits.keys():
            for neighborhood in self._results.keys():

                h_start = int(hits[hit]["Position"])
                if h_start >= self._results[neighborhood]["n_start"] and h_start <= self._results[neighborhood]["n_stop"]:
                    self._results[neighborhood]["hits"][hit] = hits[hit]
        
    def _results_filter(self, hits):
        """Add hits to the results tracker (grouped by neighbohood).
        Once all hits are added, remove neighborhood that were not
        updated.
        """
        for hit in hits.keys():
            for neighborhood in self._results.keys():

                # note that the begining and end of the hit is denoted by the begining and
                # end of the orf used as the query
                h_start = min(int(hits[hit]["Start"]), int(hits[hit]["Stop"]))
                h_stop = max(int(hits[hit]["Start"]), int(hits[hit]["Stop"]))

                # check whether this hit is contained within the neighborhood
                if h_start >= self._results[neighborhood]["n_start"] and h_stop <= self._results[neighborhood]["n_stop"]:
                    self._results[neighborhood]["hits"][hit] = hits[hit]

                    # This neighborhood contains at least one hit, so mark it
                    self._results[neighborhood]["do_not_filter"] = True 
        
        # remove unmarked neighborhood from results and query library
        remove = [neighborhood for neighborhood in self._results if not self._results[neighborhood]["do_not_filter"]]
        for neighborhood in remove:
            del self._results[neighborhood]
            del self._neighborhood_orfs[neighborhood]
        
        for neighborhood in self._results.keys():
            self._results[neighborhood]["do_not_filter"] = False

    def _get_all_orfs(self):
        """Get all of the (translated) open reading frames in this genome."""

        self._all_orfs = os.path.join(self._working_dir.name, "all_orfs.fasta")
        orffinder(sequence=self.genome, output=self._all_orfs, min_prot_len=self.min_prot_len, description=self.id)

    def _get_orfs_in_neighborhood(self, ranges):
        """Grab all of the open reading frames within a subsequence
        (neighborhood) from the original parent.   
        """
        neighborhood_orffinder(sequence=self.genome, ranges=ranges, outdir=self._working_dir.name, min_prot_len=self.min_prot_len, description=self.id)

        for r in ranges:
            key = "{}_{}".format(r[0], r[1])
            path = os.path.join(self._working_dir.name, "orf_{}_{}.fasta".format(r[0], r[1]))
            self._neighborhood_orfs[key] = path

    def add_seed_step(self, db, name, e_val, blast_type):
        """Add a seed step to the pipeline. 

        Internally, this queues a series of sub-steps that
        serve to identify genomic "neighborhoods" around the proteins of 
        interest - "seeds" - in the parent genome. 

        The individual steps can be summarized as follows:
        1. Translate all open reading frames contained in the parent genome.
        2. Blast each ORF against a reference database of target proteins (seeds)
        3. Record the span to the left and right of each hit (set by the span 
        parameter during pipeline initialization). This region is a potential 
        CRISPR-transposon gene neighborhood. Hits with overlapping regions 
        are merged into a single neighborhood.
        4. Add nighborhood and hit info to the results tracker.
        5. Remove ORFs that are not contained within the neighborhood regions
        from the query database.

        Args:
            db (str): Path to the target (seed) protein database.
            e_val (float): Blast expect value to use as a threshhold. 
            See NCBI BLAST documentation for details.
            blast_type (str): Specifies which blast program to use. 
            Currently only blastp and psiblast are supported. 

        Notes:
            Only one seed step should be added to the pipeline, and it should
            be first. Additional steps can occur in any order.
        """
        if blast_type in BLASTP_KEYWORDS:
            self._steps.append(SeedBlastp(db, name, e_val, self._working_dir.name, self.span))
        elif blast_type in PSIBLAST_KEYWORDS:
            self._steps.append(SeedBlastpsi(db, name, e_val, self._working_dir.name, self.span))
        else:
            raise ValueError("blast type option '{}' not recognized".format(blast_type))

    def add_filter_step(self, db, name, e_val, blast_type):
        """Add a filter step to the pipeline.

        Blast genomic neighborhoods against the target database. 
        Neighborhoods that don't contain hits are filtered out of 
        the results and will not be used in subsequent searches.

        Args:
            db (str): Path to the target protein database.
            e_val (float): Blast expect value to use as a threshhold. 
            See NCBI BLAST documentation for details.
            blast_type (str): Specifies which blast program to use. 
            Currently only blastp and psiblast are supported. 
        """
        if blast_type in BLASTP_KEYWORDS:
            self._steps.append(FilterBlastp(db, name, e_val, self._working_dir.name))
        elif blast_type in PSIBLAST_KEYWORDS:
            self._steps.append(FilterBlastpsi(db, name, e_val, self._working_dir.name))
        else:
            raise ValueError("blast type option '{}' not recognized".format(blast_type))
    
    def add_blast_step(self, db, name, e_val, blast_type):
        """Add a non-filtering blast step to the pipeline.

        Blast genomic neighborhoods against the target database. 
        Any hits are appended to the results.

        Args:
            db (str): Path to the target protein database.
            e_val (float): Blast expect value to use as a threshhold. 
            See NCBI BLAST documentation for details.
            blast_type (str): Specifies which blast program to use. 
            Currently only blastp and psiblast are supported.     
        """
        if blast_type in BLASTP_KEYWORDS:
            self._steps.append(Blastp(db, name, e_val, self._working_dir.name))
        elif blast_type in PSIBLAST_KEYWORDS:
            self._steps.append(Blastpsi(db, name, e_val, self._working_dir.name))
        else:
            raise ValueError("blast type option '{}' not available for filter step".format(blast_type))
    
    def add_crispr_step(self):
        """Add a pilercr step to search for CRISPR arrays."""

        self._steps.append(CrisprStep(self._working_dir.name))

    def run(self):
        """Sequentially execute each step in the pipeline.

        Currently, results from the run are returned as a dictionary
        which can, for example, be parsed or pretty-printed using 
        json:

        >>> print(json.dumps(results, indent=4))
        """

        # TODO: Add error handling for seed execution
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
                if step.is_crispr:
                    step.execute(self.genome)
                    self._results_update_crispr(step.hits)
                else:
                    step.execute(neighborhood_orfs)
                    self._results_update(step.hits)

        results = self._results
        return results


if __name__ == "__main__":

    import json

    #genome = "/home/alexis/Projects/CRISPR-Transposons/data/genomes/v_crass_J520_whole.fasta"
    genome = "/home/alexis/Projects/CRISPR-Transposons/data/contigs/C2558"
    seed_db = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/tnsAB/blast_db"
    filter_db = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/cas_uniref/blast_db"
    final_db = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/tns_dc/blast_db"

    """in_dir = "/home/alexis/Projects/CRISPR-Transposons/data/protein_references/tns_cd/"
    db_dir = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/tns_dc"
    build_blastp_db(input=in_dir, db_dir=db_dir)"""
    
    p = Pipeline(genome, "v_crass", min_prot_len=30, span=10000)
    p.add_seed_step(seed_db, "tnsAB", 0.001, "PSI")
    p.add_filter_step(filter_db, "cas", 0.001, "PROT")
    p.add_blast_step(final_db, "tnsCD", 0.001, "PROT")
    p.add_crispr_step()
    results = p.run()

    print(json.dumps(results, indent=4))
