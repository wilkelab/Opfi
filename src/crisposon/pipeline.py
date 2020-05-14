from crisposon.orffinder import orffinder, neighborhood_orffinder
from crisposon.utils import concatenate
from crisposon.build_blast_db import build_blast_db
from crisposon.steps import (SearchStep, 
                                FilterStep, 
                                SeedStep, 
                                CrisprStep, 
                                Blastp, 
                                Blastpsi, 
                                MMseqs, 
                                Diamond,
                                Pilercr)

from crisposon.output_writers import CSVWriter

import tempfile, os, json
from Bio import SeqIO

BLASTP_KEYWORDS = ["blastp", "PROT"]
PSIBLAST_KEYWORDS = ["psiblast", "PSI"]
MMSEQS_KEYWORDS = ["mmseqs"]
DIAMOND_KEYWORDS = ["diamond"]

class Pipeline:
    """
    Main class for running the CRISPR-transposon identification pipeline.
    
    Takes a single genome of interest as input, and performs a series of
    user-specified alignment steps to identify putative CRISPR-transposon 
    elements.

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

        >>> p = Pipeline(genome="v_crass.fasta", id="v_crass", span=15000)

        Add alignment steps to the pipeline. First, we will do a blast
        for TnsA/B genes - "seeds" - in the query genome. Regions that fall
        outside of the span around each hit will be filtered out of 
        subsequent searches.

        >>> p.add_blast_seed_step(db="blast_databases/tnsAB", name="tnsAB", 
                                    e_val=0.001, blast_type="PROT")

        Now add a filter step for cas proteins. This will tell the program to
        to run a blast against a cas reference database, but only for the
        regions around the hits from the previous step. 
        Regions ("neighborhoods") that do not contain hits for cas proteins
        will be filtered out of subsequent searches, as well as the results.

        >>> p.add_blast_filter_step(db="blast_databases/cas", name="cas", 
                                    e_val=0.001, blast_type="PSI")

        Finally, add a blast step for tnsC and tnsD proteins using any 
        remaining neighborhoods as queries. Neighborhoods that do not 
        contain hits will NOT be filtered out during this step.

        >>> p.add_blast_step(db="blast_databases/tnsCD", name="tnsCD", 
                                e_val=0.001, blast_type="PROT")

        Now, run the pipeline. Results are returned as a dictionary 
        object containing the hits associated with each neighborhood.

        >>> results = p.run()
    """

    def __init__(self):
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
        self._steps = []

    
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
            key = "Loc_{}-{}".format(r[0], r[1])
            self._working_results[key] = {"Loc_start-pos": int(r[0]), "Loc_end-pos": int(r[1]), 
                                    "new_hit_count": 0, "Hits": {}}

    
    def _filter(self, min_prot_count):
        """
        Remove neighborhood from results if the number of new hits added 
        during the filter step is less than the user specified cutoff.

        Note that if no cutoff is given, the default used is 1.
        """

        remove = [neighborhood for neighborhood in self._working_results 
                    if self._working_results[neighborhood]["new_hit_count"] < min_prot_count]
        
        for neighborhood in remove:
            del self._working_results[neighborhood]
            del self._neighborhood_orfs[neighborhood]
    
    
    def _results_update(self, hits, min_prot_count):
        """Add hits to the results tracker (grouped by neighborhood)."""

        for hit in hits:
            for neighborhood in self._working_results:

                # note that the begining and end of the hit is denoted by 
                # the begining and end of the orf used as the query
                h_start = min(int(hits[hit]["Query_start-pos"]), 
                                int(hits[hit]["Query_end-pos"]))
                h_stop = max(int(hits[hit]["Query_start-pos"]), 
                                int(hits[hit]["Query_end-pos"]))
                
                # check whether this hit is contained within the neighborhood
                if (h_start >= self._working_results[neighborhood]["Loc_start-pos"] 
                    and h_stop <= self._working_results[neighborhood]["Loc_end-pos"]):
                    self._working_results[neighborhood]["Hits"][hit] = hits[hit]
                    self._working_results[neighborhood]["new_hit_count"] += 1
        
        if min_prot_count >= 1:
            self._filter(min_prot_count)
        
        for neighborhood in self._working_results:
            self._working_results[neighborhood]["new_hit_count"] = 0

    
    def _results_update_crispr(self, hits):
        """Add hits for CRISPR arrays to the results tracker."""

        for hit in hits:
            for neighborhood in self._working_results:

                h_start = int(hits[hit]["Position"])
                if (h_start >= self._working_results[neighborhood]["Loc_start-pos"] 
                    and h_start <= self._working_results[neighborhood]["Loc_end-pos"]):
                    self._working_results[neighborhood]["Hits"][hit] = hits[hit]

    
    def _get_all_orfs(self, data, id):
        """Get all of the (translated) open reading frames in this genome."""

        self._all_orfs = os.path.join(self._working_dir.name, "all_orfs.fasta")
        orffinder(sequence=data, output=self._all_orfs, 
                    min_prot_len=self.min_prot_len, description=id)

    
    def _get_orfs_in_neighborhood(self, ranges, data, id):
        """
        Grab all of the open reading frames within a subsequence
        (neighborhood) from the original parent.   
        """
        self._neighborhood_orfs = {}
        neighborhood_orffinder(sequence=data, ranges=ranges, 
                                outdir=self._working_dir.name, 
                                min_prot_len=self.min_prot_len, 
                                description=id)

        for r in ranges:
            key = "Loc_{}-{}".format(r[0], r[1])
            path = os.path.join(self._working_dir.name, 
                                "orf_{}_{}.fasta".format(r[0], r[1]))
            self._neighborhood_orfs[key] = path

    
    def add_seed_step(self, db, name, e_val, blast_type, sensitivity=None, extra_args=None):
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
            self._steps.append(SeedStep(Blastp(db, e_val, name)))
        elif blast_type in PSIBLAST_KEYWORDS:
            self._steps.append(SeedStep(Blastpsi(db, e_val, name)))
        elif blast_type in MMSEQS_KEYWORDS:
            self._steps.append(SeedStep(MMseqs(db, str(e_val), name, 
                                                str(sensitivity), extra_args)))
        elif blast_type in DIAMOND_KEYWORDS:
            self._steps.append(SeedStep(Diamond(db, str(e_val), name, 
                                                str(sensitivity), extra_args)))
        else:
            raise ValueError("blast type option '{}' not recognized".format(blast_type))

    
    def add_filter_step(self, db, name, e_val, blast_type, min_prot_count=1, 
                        sensitivity=None, extra_args=None):
        """Add a filter step to the pipeline.

        Blast genomic neighborhoods against the target database. 
        Neighborhoods with no hits against the target database will 
        be filtered out of the results and will not be used in subsequent 
        searches.

        Args:
            db (str): Path to the target protein database.
            e_val (float): Blast expect value to use as a threshhold. 
                See NCBI BLAST documentation for details.
            blast_type (str): Specifies which blast program to use. 
                Currently only blastp and psiblast are supported.
            min_prot_count (int, optional): Sets a minimum number of
                hits required to keep neighborhoods. 
        """
        if blast_type in BLASTP_KEYWORDS:
            self._steps.append(FilterStep(Blastp(db, e_val, name), min_prot_count))
        elif blast_type in PSIBLAST_KEYWORDS:
            self._steps.append(FilterStep(Blastpsi(db, e_val, name), min_prot_count))
        elif blast_type in MMSEQS_KEYWORDS:
            self._steps.append(FilterStep(MMseqs(db, str(e_val), name, str(sensitivity), 
                                                    extra_args), min_prot_count))
        elif blast_type in DIAMOND_KEYWORDS:
            self._steps.append(FilterStep(Diamond(db, str(e_val), name, str(sensitivity), 
                                                    extra_args), min_prot_count))
        else:
            raise ValueError("blast type option '{}' not recognized".format(blast_type))
    
    
    def add_blast_step(self, db, name, e_val, blast_type, 
                        sensitivity=None, extra_args=None):
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
            self._steps.append(SearchStep(Blastp(db, e_val, name)))
        elif blast_type in PSIBLAST_KEYWORDS:
            self._steps.append(SearchStep(Blastpsi(db, e_val, name)))
        elif blast_type in MMSEQS_KEYWORDS:
            self._steps.append(SearchStep(MMseqs(db, str(e_val), name, str(sensitivity), extra_args)))
        elif blast_type in DIAMOND_KEYWORDS:
            self._steps.append(SearchStep(Diamond(db, str(e_val), name, str(sensitivity), extra_args)))
        else:
            raise ValueError("blast type option '{}' not available for filter step".format(blast_type))
    
    
    def add_crispr_step(self):
        """Add a step to search for CRISPR arrays.
        
        Uses pilercr with default parameters. Hits that 
        overlap with a genomic neighborhood are appended to
        the resutls.
        """

        self._steps.append(CrisprStep(Pilercr("CRISPR")))

    
    def _format_results(self, outfrmt, outfile):
        
        # Remove temporary hit counter tag
        for contig in self._results:
            for neighborhood in self._results[contig]:
                del self._results[contig][neighborhood]["new_hit_count"]
        
        if outfrmt is not None:
            if outfrmt == "JSON":
                try:
                    with open(outfile, 'w') as jsonfile:
                        json.dump(self._results, jsonfile)
                except FileNotFoundError:
                    with open("results.json", 'w') as jsonfile:
                        json.dump(self._results, jsonfile)

            elif outfrmt == "CSV":
                csv_writer = CSVWriter(self._results, outfile)
                csv_writer.to_csv()
    
    
    def _record_all_hits(self, outfile):
        """Write intermediate hits to a json file."""
        
        try:
            with open(outfile, "w") as jsonfile:
                json.dump(self._all_hits, jsonfile)
        
        except (FileNotFoundError, TypeError) as e:
            with open("all_hits.json", "w") as jsonfile:
                json.dump(self._all_hits, jsonfile)
            
            if isinstance(e, FileNotFoundError):
                print("Cannot open {}".format(outfile),
                        " writing all hits to working directory")
            else:
                print("No output file given for writing all" +
                        " hits, using working directory")

    
    def run(self, data, min_prot_len=30, span=10000,
            outfrmt=None, outfile=None, record_all_hits=False,
            all_hits_outfile=None):
        """Sequentially execute each step in the pipeline.

        Currently, results from the run are returned as a dictionary
        which can, for example, be parsed or pretty-printed using 
        json:

        >>> print(json.dumps(results, indent=4))
        """
        
        self.data_path = data
        self.min_prot_len = min_prot_len
        self.span = span

        self._results = {}
        self._all_hits = {}

        for record in SeqIO.parse(data, "fasta"):
            
            contig_id = record.id
            self._working_dir = tempfile.TemporaryDirectory()
            contig_path = os.path.join(self._working_dir.name, "contig.fasta")
            SeqIO.write(record, contig_path, "fasta")

            self._get_all_orfs(contig_path, contig_id)
            self._working_results = {}
            self._all_hits[contig_id] = {}

            neighborhood_orfs = None
            for step in self._steps:
                
                if isinstance(step, SeedStep):
                    #print("Begin seed step: {}".format(step.name))
                    step.execute(self._all_orfs, self.span)

                    if len(step.hits) != 0:
                        self._get_orfs_in_neighborhood(step.neighborhood_ranges, 
                                                        contig_path,
                                                        contig_id)
                        self._results_init(step.neighborhood_ranges)
                        self._results_update(step.hits, min_prot_count=0)
                        neighborhood_orfs = concatenate(self._working_dir.name, 
                                                        self._neighborhood_orfs.values())
                        #print("Seed step complete")
                    
                    else:
                        #print("No hits for seed gene - terminating run")
                        self._results[contig_id] = {}
                        break

                elif isinstance(step, FilterStep):
                    #print("Begin filter step: {}".format(step.name))
                    
                    if len(self._neighborhood_orfs) != 0:
                        step.execute(neighborhood_orfs)
                        self._results_update(step.hits, min_prot_count=step.min_prot_count)
                        neighborhood_orfs = concatenate(self._working_dir.name, 
                                                        self._neighborhood_orfs.values())
                        #print("Filtering for {} genes complete".format(step.name))
                    
                    else:
                        #print("No putative neighborhoods remain - terminating run")
                        self._results[contig_id] = {}
                        break
                
                elif isinstance(step, CrisprStep):
                    #print("Begin CRISPR array search")
                    
                    if len(self._neighborhood_orfs) != 0:
                        step.execute(contig_path)
                        self._results_update_crispr(step.hits)
                        #print("CRISPR array search complete")
                    
                    else:
                        #print("No putative neighborhoods remain - terminating run")
                        self._results[contig_id] = {}
                        break
                        
                else:
                    #print("Begin blast step: {}".format(step.name))
                    if len(self._neighborhood_orfs) != 0:
                        step.execute(neighborhood_orfs)
                        self._results_update(step.hits, min_prot_count=0)
                        #print("Blast step complete")
                    
                    else:
                        #print("No putative neighborhoods remain - terminating run")
                        self._results[contig_id] = {}
                        break
                
                if record_all_hits:
                    self._all_hits[contig_id][step.search_tool.step_id] = step.hits
                
                self._results[contig_id] = self._working_results
        
        self._format_results(outfrmt, outfile)
        results = self._results

        if record_all_hits:
            self._record_all_hits(all_hits_outfile)
        
        return results
