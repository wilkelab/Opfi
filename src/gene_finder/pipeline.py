from gene_finder.orffinder import orffinder, neighborhood_orffinder
from gene_finder.utils import concatenate
from gene_finder.steps import (SearchStep, 
                                FilterStep, 
                                SeedStep,
                                SeedWithCoordinatesStep, 
                                CrisprStep, 
                                Blastp, 
                                Blastpsi, 
                                MMseqs, 
                                Diamond,
                                Pilercr)

from gene_finder.output_writers import CSVWriter

import tempfile, os, json, gzip
from Bio import SeqIO

class Pipeline:
    """
    Main class for running the CRISPR-transposon identification pipeline.
    
    Performs a series of user-specified local alignment steps to identify 
    putative CRISPR-transposon elements.

    Example:

        Create a pipeline to search for CRISPR-transposon systems
        in the Vibrio crassostreae genome.

        >>> p = Pipeline()

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

        >>> results = p.run(data="v_crassostreae.fasta")
    """

    def __init__(self):
        """Initialize a Pipeline object."""
        # defined in pipeline.run()
        self.data_path = None
        self.min_prot_len = None
        self.span = None
        
        # modified by add_step methods
        self._steps = []
        
        # collect results about the input data
        # will be reset each time pipeline.run() is called
        self._results = {}
        self._all_hits = {}

        # collect data specific to a single contig in the input
        # will be reset each time a new contig is processed
        self._working_results = {}
        self._all_orfs = None
        self._neighborhood_orfs = {}
        self._working_dir = None

    
    def __del__(self):
        """
        Delete the working directory and its contents when 
        this object is garbage collected.
        """
        if self._working_dir is not None:
            self._working_dir.cleanup()
    
    def _setup_run(self, record):
        """
        Prepare for running the pipeline on a new contig.
        Since the same setup is re-used on multple contigs
        (or even multiple fasta files) we want to make sure that
        old data is cleared out.

        Also gather info about the current contig and run 
        orffinder to get all ORFs in this contig.
        """
        self._reset_contig_data()
        self._working_dir = tempfile.TemporaryDirectory()
        contig_id = record.id
        contig_len = len(record.seq)
        contig_path = os.path.join(self._working_dir.name, "contig.fasta")
        SeqIO.write(record, contig_path, "fasta")
        self._get_all_orfs(contig_path, contig_id)
        self._all_hits[contig_id] = {}
        return contig_id, contig_len, contig_path
    
    def _reset_contig_data(self):
        """Clear out data from a previous contig."""
        if self._working_dir is not None:
            self._working_dir.cleanup()
        self._working_results = {}
        self._neighborhood_orfs = {}
        self._all_orfs = None
    
    def _reset_results(self):
        """Clear out results from a previous run."""
        self._results = {}
        self._all_hits = {}

    def _open_data(self, data, is_binary):
        """
        Return the appropriate handle/file object for reading
        input data that is either flat text or a gzipped file. 
        """
        if is_binary:
            handle = gzip.open(data, "rt")
        else:
            handle = open(data, "r")
        
        return handle

    def _results_init(self, neighborhood_ranges):
        """
        Create an entry for every gene neighborhood identified
        during the seed phase in the working results tracker.
        """
        for r in neighborhood_ranges:
            key = "Loc_{}-{}".format(r[0], r[1])
            self._working_results[key] = {"Loc_start-pos": int(r[0]), "Loc_end-pos": int(r[1]), 
                                    "new_hit_count": 0, "Hits": {}}

    
    def _filter(self, min_prot_count):
        """
        Remove neighborhood from the working results if the number of
        new hits added during the filter step is less than the user 
        specified cutoff.

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

        orfs = os.path.join(self._working_dir.name, "all_orfs.fasta")
        self._all_orfs = orffinder(sequence=data, output=orfs, 
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

    
    def add_seed_step(self, db, name, e_val, blast_type, sensitivity=None, parse_descriptions=True, **kwargs):
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
            name (str): A unique name/ID for this step in the pipeline.
            e_val (float): Expect value to use as a threshhold. 
            blast_type (str): Specifies which search program to use. 
                This can be either "PROT" (blastp), "PSI" (psiblast),
                "mmseqs" (mmseqs2), or "diamond" (diamond). Note that
                mmseqs2 and diamond support is currently experimental.
            sensitivity (str): Sets the sensitivity param 
                for mmseqs and diamond (does nothing if blast is the
                seach type).
            parse_descriptions (bool, optional): By default, reference protein
                descriptions (from fasta headers) are parsed for gene name labels;
                specifically, descriptions are split on whitespace characters
                and the second item is used for the label. Make this false to 
                simply use the whole protein description for the label 
                (i.e everything after the first whitespace in the header). If
                using this option with NCBI blast, DO NOT use the `-parse_seqids`
                flag when creating protein databases with `makeblastdb`.
            **kwargs: These can be any additional blast parameters,
                specified as key-value pairs. Note that certain parameters
                are not allowed, mainly those that control output formatting.
                Currently only supported for blastp/psiblast; if blast_type
                is set to mmseqs or diamond, kwargs will be silently ignored.

        Notes:
            Only one seed step should be added to the pipeline, and it should
            be first. Additional steps can occur in any order.
        """
        if blast_type == "PROT" or "blastp":
            self._steps.append(SeedStep(Blastp(db, e_val, name, parse_descriptions, kwargs)))
        elif blast_type == "PSI" or "psiblast":
            self._steps.append(SeedStep(Blastpsi(db, e_val, name, parse_descriptions, kwargs)))
        elif blast_type == "mmseqs":
            self._steps.append(SeedStep(MMseqs(db, str(e_val), name, str(sensitivity), parse_descriptions)))
        elif blast_type == "diamond":
            self._steps.append(SeedStep(Diamond(db, str(e_val), name, str(sensitivity), parse_descriptions)))
        else:
            raise ValueError("blast type option '{}' not recognized".format(blast_type))
    

    def add_seed_with_coordinates_step(self, db, e_val, blast_type, title=None, sensitivity=None, 
                                       start=None, end=None, contig_id=None, **kwargs):
        """
        Define a genomic region to search with coordinates instead of a bait gene.

        An alternative to `add_seed_step`. Most useful for re-annotating 
        putative systems of interest, where the region coordinates are already
        known.

        Args:
            db (str): Path to the target database to search against.
            e_val (float): Expect value to use as a threshhold. 
            blast_type (str): Specifies which search program to use. 
                This can be either "PROT" (blastp), "PSI" (psiblast),
                "mmseqs" (mmseqs2), or "diamond" (diamond). Note that
                mmseqs2 and diamond support is currently experimental.
            title (str, optional): A unique name for this step. If no 
                title is given, the pipeline automatically assigns one
                based on the order the step was added, i.e. "step_<index>"
            sensitivity (str): Sets the sensitivity param 
                for mmseqs and diamond (does nothing if blast is the
                seach type).
            start (int): Defines the beginning of the region to search (in bp).
                If no start position is given the first (zero indexed) position
                in the contig is used.
            end (int): Defines the end of the region to search (in bp). If no
                end position is given the last position in the contig is used.
            contig_id(string, optional): An identifier for the contig to search.
                If no ID is given, the pipeline will search every contig in the
                input file using the coordinates specified. Note that the contig ID 
                is defined as the substring between the ">" character and the first
                " " character in the contig header.
        """
        self._steps.append(SeedWithCoordinatesStep(start=start, end=end, contig_id=contig_id))

    def add_filter_step(self, db, name, e_val, blast_type, min_prot_count=1, 
                        sensitivity=None, parse_descriptions=True, **kwargs):
        """Add a filter step to the pipeline.

        Blast genomic neighborhoods against the target database. 
        Neighborhoods with no hits against the target database will 
        be filtered out of the results and will not be used in subsequent 
        searches.

        Args:
            db (str): Path to the target (seed) protein database.
            name (str): A unique name/ID for this step in the pipeline.
            e_val (float): Expect value to use as a threshhold. 
            blast_type (str): Specifies which search program to use. 
                This can be either "PROT" (blastp), "PSI" (psiblast),
                "mmseqs" (mmseqs2), or "diamond" (diamond). Note that
                mmseqs2 and diamond support is currently experimental.
            min_prot_count (int, optional): Minimum number of hits 
                needed for the neighborhood to be retained. Default 
                is one.
            sensitivity (str): Sets the sensitivity param 
                for mmseqs and diamond (does nothing if blast is the
                seach type).
            parse_descriptions (bool, optional): By default, reference protein
                descriptions (from fasta headers) are parsed for gene name labels;
                specifically, descriptions are split on whitespace characters
                and the second item is used for the label. Make this false to 
                simply use the whole protein description for the label 
                (i.e everything after the first whitespace in the header). If
                using this option with NCBI blast, DO NOT use the `-parse_seqids`
                flag when creating protein databases with `makeblastdb`.  
            **kwargs: These can be any additional blast parameters,
                specified as key-value pairs. Note that certain parameters
                are not allowed, mainly those that control output formatting.
                Currently only supported for blastp/psiblast; if blast_type
                is set to mmseqs or diamond, kwargs will be silently ignored.
        """
        if blast_type == "PROT" or "blastp":
            self._steps.append(FilterStep(Blastp(db, e_val, name, parse_descriptions, kwargs), min_prot_count))
        elif blast_type == "PSI" or "psiblast":
            self._steps.append(FilterStep(Blastpsi(db, e_val, name, parse_descriptions, kwargs), min_prot_count))
        elif blast_type == "mmseqs":
            self._steps.append(FilterStep(MMseqs(db, str(e_val), name, str(sensitivity), parse_descriptions), min_prot_count))
        elif blast_type == "diamond":
            self._steps.append(FilterStep(Diamond(db, str(e_val), name, str(sensitivity), parse_descriptions), min_prot_count))
        else:
            raise ValueError("blast type option '{}' not recognized".format(blast_type))
    
    
    def add_blast_step(self, db, name, e_val, blast_type, 
                        sensitivity=None, parse_descriptions=True, **kwargs):
        """Add a non-filtering blast step to the pipeline.

        Blast genomic neighborhoods against the target database. 
        Any hits are appended to the results.

        Args:
            db (str): Path to the target (seed) protein database.
            name (str): A unique name/ID for this step in the pipeline.
            e_val (float): Expect value to use as a threshhold. 
            blast_type (str): Specifies which search program to use. 
                This can be either "PROT" (blastp), "PSI" (psiblast),
                "mmseqs" (mmseqs2), or "diamond" (diamond). Note that
                mmseqs2 and diamond support is currently experimental.
            min_prot_count (int, optional): Minimum number of hits 
                needed for the neighborhood to be retained. Default 
                is one.
            sensitivity (str): Sets the sensitivity param 
                for mmseqs and diamond (does nothing if blast is the
                seach type). 
            parse_descriptions (bool, optional): By default, reference protein
                descriptions (from fasta headers) are parsed for gene name labels;
                specifically, descriptions are split on whitespace characters
                and the second item is used for the label. Make this false to 
                simply use the whole protein description for the label 
                (i.e everything after the first whitespace in the header). If
                using this option with NCBI blast, DO NOT use the `-parse_seqids`
                flag when creating protein databases with `makeblastdb`.
            **kwargs: These can be any additional blast parameters,
                specified as key-value pairs. Note that certain parameters
                are not allowed, mainly those that control output formatting.
                Currently only supported for blastp/psiblast; if blast_type
                is set to mmseqs or diamond, kwargs will be silently ignored.     
        """
        if blast_type == "PROT" or "blastp":
            self._steps.append(SearchStep(Blastp(db, e_val, name, parse_descriptions, kwargs)))
        elif blast_type == "PSI" or "psiblast":
            self._steps.append(SearchStep(Blastpsi(db, e_val, name, parse_descriptions, kwargs)))
        elif blast_type == "mmseqs":
            self._steps.append(SearchStep(MMseqs(db, str(e_val), name, str(sensitivity), parse_descriptions)))
        elif blast_type == "diamond":
            self._steps.append(SearchStep(Diamond(db, str(e_val), name, str(sensitivity), parse_descriptions)))
        else:
            raise ValueError("blast type option '{}' not available for filter step".format(blast_type))
    
    
    def add_crispr_step(self):
        """Add a step to search for CRISPR arrays.
        
        Uses pilercr with default parameters. Hits that 
        overlap with a genomic neighborhood are appended to
        the resutls.
        """

        self._steps.append(CrisprStep(Pilercr("CRISPR")))
    

    def _update_output_sequences(self):
        for neighborhood, path in self._neighborhood_orfs.items():
            sequences = {}
            for record in SeqIO.parse(path, "fasta"):
                sequences[record.id] = str(record.seq)
            for key, hit in self._working_results[neighborhood]["Hits"].items():
                if "Query_ORFID" in hit.keys():
                    self._working_results[neighborhood]["Hits"][key]["Query_seq"] = sequences[hit["Query_ORFID"]]

    
    def _format_results(self, outfrmt, outfile):
        """Process results into their final format.

        If an output format was specified, also 
        writes results to either a JSON file or
        a CSV file.
        """
        
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
                csv_writer.to_csv(self.data_path)
    
    
    def _record_all_hits(self, outfile):
        """Write all/intermediate hits to a json file."""
        
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

    
    def run(self, data, min_prot_len=60, span=10000,
            outfrmt=None, outfile=None, record_all_hits=False,
            all_hits_outfile=None, gzip=False) -> dict:
        """Sequentially execute each step in the pipeline.

        Args:
            data (str): Path to input data file. Can be a single-
                or multi-sequence file in fasta format.
            min_prot_len (int, optional): Minimum ORF length (aa).
                Default is 60.
            span (int, optional): Length (nt) upsteam and downstream
                of each seed hit to keep. Defines the aproximate size
                of the genomic neighborhoods that will be used as the
                search space after the seed step.
            outfrmt (str, optional): Specifies the output file format.
                Can be either "CSV" or "JSON". If no output format is
                given then the results will not be written to disk.
            outfile (str, optional): Path to the file to write results
                to.
            record_all_hits (bool, optional): If set to True then 
                all hits against all references will be written to
                a file.
            all_hits_outfile (str, optional): Path to the file to
                write all hit data to.

        Returns:   
            Results (dict): Candidate systems, grouped by contig id
                and genomic location.
        """
        self._reset_results()
        self.data_path = data
        self.min_prot_len = min_prot_len
        self.span = span

        data_handle = self._open_data(self.data_path, gzip)
        for record in SeqIO.parse(data_handle, "fasta"):
            # if this run was seeded with a particular contig id, skip all other contigs
            # in the file if they exit
            # saves a bit of time, since there is overhead associated with getting things set up
            if isinstance(self._steps[0], SeedWithCoordinatesStep):
                if self._steps[0].contig_id is not None and self._steps[0].contig_id != record.id:
                    continue
            # clear out any data that may be leftover from processing a previous
            # contig and get all ORFs in this contig
            contig_id, contig_len, contig_path = self._setup_run(record)
            if self._all_orfs is None:
                # No ORFs were identified in the contig (probably because it's small),
                # so we just continue on to the next contig if it exists
                self._all_hits[contig_id] = {}
                self._results[contig_id] = {}
                continue
            for step in self._steps:
                if isinstance(step, SeedStep):
                    #print("Begin seed step: {}".format(step.name))
                    step.execute(self._all_orfs, self.span, contig_len)
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

                elif isinstance(step, SeedWithCoordinatesStep):
                    if step.start is None:
                        step.update_start_coord(0)
                    if step.end is None:
                        step.update_end_coord(contig_len)
                    self._get_orfs_in_neighborhood(step.neighborhood_ranges,
                                                   contig_path,
                                                   contig_id)
                    self._results_init(step.neighborhood_ranges)
                    neighborhood_orfs = concatenate(self._working_dir.name,
                                                    self._neighborhood_orfs.values())
                    # skip hit processing at the end of the loop, since we don't have hits yet
                    continue

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

                self._all_hits[contig_id][step.search_tool.step_id] = step.hits
                self._update_output_sequences()
                self._results[contig_id] = self._working_results
        
        self._format_results(outfrmt, outfile)
        if record_all_hits:
            self._record_all_hits(all_hits_outfile)
        
        data_handle.close() # make sure to close the input data file object
        return self._results
