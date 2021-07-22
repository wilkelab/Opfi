from gene_finder.orffinder import orffinder, neighborhood_orffinder
from gene_finder.utils import concatenate
from gene_finder.steps import (SearchStep, 
                                FilterStep, 
                                SeedStep,
                                SeedWithCoordinatesStep, 
                                CrisprStep, 
                                Blastn,
                                BlastnStep,
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
    Coordinates protein (or nucleic acid) searches to find gene clusters
    of interest in genomic/metagenomic data.
    """

    def __init__(self):
        """
        Create a pipeline for finding gene cluster candidates.
        """
        # Most members are re-assigned later, either when search steps are added 
        # or when Pipeline.run() is called, but each is intialized here for clarity

        # defined in Pipeline.run()
        self.data_path = None
        self.min_prot_len = None
        self.span = None
        self.job_id = None
        self.output_directory = None
        
        # modified by add_step methods
        self._steps = []
        
        # collect results about the input data
        # will be reset each time Pipeline.run() is called
        self._results = {}
        self._all_hits = {}

        # basically give the current state of the output files; that is, 
        # whether any data from this run has been written to disk yet
        # these are also reset (to this state) when Pipeline.run() is called
        self._appending_results = False
        self._appending_hits = False

        # collect data specific to a single contig in the input
        # will be reset each time a new contig is processed
        self._working_results = {}
        self._all_orfs = None
        self._neighborhood_orfs = {}
        self._working_dir = None

    
    def __del__(self):
        """
        Delete the working directory and its contents when this object is garbage collected.
        """
        if self._working_dir is not None:
            self._working_dir.cleanup()
    

    def _setup_run(self, record):
        """
        Prepare for running the pipeline on a new contig. Since the same setup is re-used on 
        multple contigs (or even multiple fasta files) we want to make sure that old data is 
        cleared out.

        Also gather info about the current contig and call `orffinder` to get all ORFs in
        this contig.
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
        """
        Clear out data from a previous contig.
        """
        if self._working_dir is not None:
            self._working_dir.cleanup()
        self._working_results = {}
        self._neighborhood_orfs = {}
        self._all_orfs = None
    

    def _reset_results(self):
        """
        Clear out results from a previous run.
        """
        self._results = {}
        self._all_hits = {}
        self._appending_hits = False
        self._appending_results = False


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
        Remove neighborhoods from the working results if the number of
        new hits added during the filter step is less than the user 
        specified cutoff. If no cutoff is given the default used is 1.
        """
        remove = [neighborhood for neighborhood in self._working_results 
                    if self._working_results[neighborhood]["new_hit_count"] < min_prot_count]
        
        for neighborhood in remove:
            del self._working_results[neighborhood]
            del self._neighborhood_orfs[neighborhood]
    
    
    def _results_update(self, hits, min_prot_count):
        """
        Add hits to the results tracker (grouped by neighborhood).
        """
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
        """
        Add hits for CRISPR arrays to the results tracker.
        """
        for hit in hits:
            for neighborhood in self._working_results:

                h_start = int(hits[hit]["Position"])
                if (h_start >= self._working_results[neighborhood]["Loc_start-pos"] 
                    and h_start <= self._working_results[neighborhood]["Loc_end-pos"]):
                    self._working_results[neighborhood]["Hits"][hit] = hits[hit]

    
    def _results_update_nucl(self, hits):
        """
        Add hits for nucleotide BLAST.
        """
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
        
        for neighborhood in self._working_results:
            self._working_results[neighborhood]["new_hit_count"] = 0


    def _get_all_orfs(self, data, id):
        """
        Get all of the (translated) open reading frames in this genome.
        """
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

    
    def add_seed_step(self, db, name, e_val, blast_type, sensitivity=None, parse_descriptions=True, blast_path=None, **kwargs):
        """
        Find genomic regions that contain at least one "seed" sequence.

        Args:
            db (str): Path to the target (seed) protein database.
            name (str): A unique name/ID for this step in the pipeline.
            e_val (float): Expect value to use. Only keep hits with a an equivalent or
                better (lower) score.
            blast_type (str): Specifies which search program to use. 
                This can be either "PROT" (blastp), "PSI" (psiblast),
                "mmseqs" (mmseqs2), or "diamond" (diamond).
            sensitivity (str): Sets the sensitivity param 
                for mmseqs and diamond (does nothing if BLAST is the
                seach type).
            parse_descriptions (bool, optional): By default, reference protein
                descriptions (from fasta headers) are parsed for gene name labels;
                specifically, descriptions are split on whitespace characters
                and the second item is used for the label. Make this false to 
                simply use the whole protein description for the label 
                (i.e everything after the first whitespace in the header). If
                using this option with NCBI BLAST, DO NOT use the ``-parse_seqids``
                flag when creating protein databases with :program:`makeblastdb`.
            blast_path (string, optional): Path to the blastp/mmseqs/diamond program,
                if not using the system default.
            **kwargs: These can be any additional BLAST parameters,
                specified as key-value pairs. Note that certain parameters
                are not allowed, mainly those that control output formatting.
                Currently only supported for blastp/psiblast; if blast_type
                is set to mmseqs or diamond, kwargs will be silently ignored.

        Note:
            This should be the first step added to a :class:`gene_finder.pipeline.Pipeline` object.
            Additional gene finding steps can be added in any order. 
        """
        if blast_type in ("PROT", "blastp"):
            path = 'blastp' if blast_path is None else blast_path
            self._steps.append(SeedStep(Blastp(db, e_val, name, parse_descriptions, path, kwargs)))
        elif blast_type in ("PSI", "psiblast"):
            path = 'psiblast' if blast_path is None else blast_path
            self._steps.append(SeedStep(Blastpsi(db, e_val, name, parse_descriptions, path, kwargs)))
        elif blast_type == "mmseqs":
            self._steps.append(SeedStep(MMseqs(db, str(e_val), name, sensitivity, parse_descriptions)))
        elif blast_type == "diamond":
            self._steps.append(SeedStep(Diamond(db, str(e_val), name, sensitivity, parse_descriptions)))
        else:
            raise ValueError("blast type option '{}' not recognized".format(blast_type))
    

    def add_seed_with_coordinates_step(self, db, name, e_val, blast_type, sensitivity=None, 
                                       parse_descriptions=True, start=None, end=None, 
                                       contig_id=None, blast_path=None, **kwargs):
        """
        Define a genomic region of interest with coordinates instead of a seed sequence.

        An alternative to :meth:`gene_finder.pipeline.Pipeline.add_seed_step`. Most useful for re-annotating 
        putative systems of interest, where the region coordinates are already
        known.

        Args:
            db (str): Path to the target database to search against.
            name (str): A unique name/ID for this step in the pipeline.
            e_val (float): Expect value to use. Only keep hits with an equivalent or
                better (lower) score.
            blast_type (str): Specifies which search program to use. 
                This can be either "PROT" (blastp), "PSI" (psiblast),
                "mmseqs" (mmseqs2), or "diamond" (diamond).
            sensitivity (str): Sets the sensitivity param 
                for mmseqs and diamond (does nothing if BLAST is the
                seach type).
            parse_descriptions (bool, optional): By default, reference protein
                descriptions (from fasta headers) are parsed for gene name labels;
                specifically, descriptions are split on whitespace characters
                and the second item is used for the label. Make this false to 
                simply use the whole protein description for the label 
                (i.e everything after the first whitespace in the header). If
                using this option with NCBI BLAST, DO NOT use the ``-parse_seqids``
                flag when creating protein databases with :program:`makeblastdb`.
            start (int): Defines the beginning of the region to search, in base pairs (bp).
                If no start position is given the first (zero indexed) position
                in the genome/contig is used.
            end (int): Defines the end of the region to search, in base pairs (bp). If no
                end position is given the last position in the contig is used.
            contig_id(string, optional): An identifier for the contig to search.
                If no ID is given, the pipeline will search every contig in the
                input file using the coordinates specified. Note that the contig ID 
                is defined as the substring between the ">" character and the first
                " " character in the contig header.
            blast_path (string, optional): Path to the blastp/mmseqs/diamond program,
                if not using the system default.
            **kwargs: These can be any additional BLAST parameters,
                specified as key-value pairs. Note that certain parameters
                are not allowed, mainly those that control output formatting.
                Currently only supported for blastp/psiblast; if blast_type
                is set to mmseqs or diamond, kwargs will be silently ignored.
        """
        self._steps.append(SeedWithCoordinatesStep(start=start, end=end, contig_id=contig_id))
        self.add_blast_step(db=db, name=name, e_val=e_val, blast_type=blast_type, 
                            sensitivity=sensitivity, parse_descriptions=parse_descriptions, blast_path=blast_path, **kwargs)


    def add_filter_step(self, db, name, e_val, blast_type, min_prot_count=1, 
                        sensitivity=None, parse_descriptions=True, blast_path=None, **kwargs):
        """
        Add a step to search candidate regions for target sequences, and filter out
        candidates that do not have at least ``min_prot_count`` matching sequences.

        Args:
            db (str): Path to the target protein sequence database.
            name (str): A unique name/ID for this step in the pipeline.
            e_val (float): Expect value to use. Only keep hits with a an equivalent or
                better (lower) score.
            blast_type (str): Specifies which search program to use. 
                This can be either "PROT" (blastp), "PSI" (psiblast),
                "mmseqs" (mmseqs2), or "diamond" (diamond).
            min_prot_count (int, optional): Minimum number of hits 
                needed to keep each candidate.
            sensitivity (str): Sets the sensitivity param 
                for mmseqs and diamond (does nothing if BLAST is the
                seach type).
            parse_descriptions (bool, optional): By default, reference protein
                descriptions (from fasta headers) are parsed for gene name labels;
                specifically, descriptions are split on whitespace characters
                and the second item is used for the label. Make this false to 
                simply use the whole protein description for the label 
                (i.e everything after the first whitespace in the header). If
                using this option with NCBI blast, DO NOT use the ``-parse_seqids``
                flag when creating protein databases with :program:`makeblastdb`.  
            blast_path (string, optional): Path to the blastp/mmseqs/diamond program,
                if not using the system default.
            **kwargs: These can be any additional BLAST parameters,
                specified as key-value pairs. Note that certain parameters
                are not allowed, mainly those that control output formatting.
                Currently only supported for blastp/psiblast; if blast_type
                is set to mmseqs or diamond, kwargs will be silently ignored.
        """
        if blast_type in ("PROT", "blastp"):
            path = 'blastp' if blast_path is None else blast_path
            self._steps.append(FilterStep(Blastp(db, e_val, name, parse_descriptions, path, kwargs), min_prot_count))
        elif blast_type in ("PSI", "psiblast"):
            path = 'psiblast' if blast_path is None else blast_path
            self._steps.append(FilterStep(Blastpsi(db, e_val, name, parse_descriptions, path, kwargs), min_prot_count))
        elif blast_type == "mmseqs":
            self._steps.append(FilterStep(MMseqs(db, str(e_val), name, sensitivity, parse_descriptions), min_prot_count))
        elif blast_type == "diamond":
            self._steps.append(FilterStep(Diamond(db, str(e_val), name, sensitivity, parse_descriptions), min_prot_count))
        else:
            raise ValueError("blast type option '{}' not recognized".format(blast_type))
    
    
    def add_blast_step(self, db, name, e_val, blast_type, 
                        sensitivity=None, parse_descriptions=True, blast_path=None, **kwargs):
        """
        Add a non-filtering search step to the pipeline. That is, search each candidate 
        for target sequences without applying any filtering logic. This is most useful
        for annotating candidates for non-essential or ancillary genes.

        Args:
            db (str): Path to the target protein sequence database.
            name (str): A unique name/ID for this step in the pipeline.
            e_val (float): Expect value to use. Only keep hits with a an equivalent or
                better (lower) score.
            blast_type (str): Specifies which search program to use. 
                This can be either "PROT" (blastp), "PSI" (psiblast),
                "mmseqs" (mmseqs2), or "diamond" (diamond).
            sensitivity (str): Sets the sensitivity param 
                for mmseqs and diamond (does nothing if BLAST is the
                seach type). 
            parse_descriptions (bool, optional): By default, reference protein
                descriptions (from fasta headers) are parsed for gene name labels;
                specifically, descriptions are split on whitespace characters
                and the second item is used for the label. Make this false to 
                simply use the whole protein description for the label 
                (i.e everything after the first whitespace in the header). If
                using this option with NCBI BLAST, DO NOT use the ``-parse_seqids``
                flag when creating protein databases with :program:`makeblastdb`.
            blast_path (string, optional): Path to the blastp/mmseqs/diamond program,
                if not using the system default.
            **kwargs: These can be any additional BLAST parameters,
                specified as key-value pairs. Note that certain parameters
                are not allowed, mainly those that control output formatting.
                Currently only supported for blastp/psiblast; if blast_type
                is set to mmseqs or diamond, kwargs will be silently ignored.     
        """
        if blast_type in ("PROT", "blastp"):
            path = 'blastp' if blast_path is None else blast_path
            self._steps.append(SearchStep(Blastp(db, e_val, name, parse_descriptions, path, kwargs)))
        elif blast_type in ("PSI", "psiblast"):
            path = 'psiblast' if blast_path is None else blast_path
            self._steps.append(SearchStep(Blastpsi(db, e_val, name, parse_descriptions, path, kwargs)))
        elif blast_type == "mmseqs":
            self._steps.append(SearchStep(MMseqs(db, str(e_val), name, sensitivity, parse_descriptions)))
        elif blast_type == "diamond":
            self._steps.append(SearchStep(Diamond(db, str(e_val), name, sensitivity, parse_descriptions)))
        else:
            raise ValueError("blast type option '{}' not available for filter step".format(blast_type))
    
    
    def add_crispr_step(self):
        """
        Add a step to search for CRISPR arrays using PILER-CR.
        """
        self._steps.append(CrisprStep(Pilercr("CRISPR")))
    

    def add_blastn_step(self, db, name, e_val, parse_descriptions=False, blastn_path='blastn', **kwargs):
        """ 
        Add a step to do nucleotide BLAST. 

        Args:
            db (str): Path to the target protein sequence database.
            name (str): A unique name/ID for this step in the pipeline.
            e_val (float): Expect value to use. Only keep hits with a an equivalent or
                better (lower) score.
            parse_descriptions (bool, optional): By default, reference protein
                descriptions (from fasta headers) are parsed for gene name labels;
                specifically, descriptions are split on whitespace characters
                and the second item is used for the label. Make this false to 
                simply use the whole protein description for the label 
                (i.e everything after the first whitespace in the header). If
                using this option with NCBI BLAST, DO NOT use the ``-parse_seqids``
                flag when creating protein databases with :program:`makeblastdb`.
            blast_path (string, optional): Path to the blastn program,
                if not using the system default.
            **kwargs: These can be any additional BLAST parameters,
                specified as key-value pairs. Note that certain parameters
                are not allowed, mainly those that control output formatting.
                Currently only supported for blastp/psiblast; if blast_type
                is set to mmseqs or diamond, kwargs will be silently ignored.
        """
        self._steps.append(BlastnStep(Blastn(db, name, e_val, parse_descriptions, blastn_path, kwargs)))


    def _update_output_sequences(self):
        """
        Replace alignment sequences (which contain gap characters) with full ORF sequences 
        in the final output.
        """
        for neighborhood, path in self._neighborhood_orfs.items():
            sequences = {}
            for record in SeqIO.parse(path, "fasta"):
                sequences[record.id] = str(record.seq)
            for key, hit in self._working_results[neighborhood]["Hits"].items():
                if "Query_ORFID" in hit.keys():
                    self._working_results[neighborhood]["Hits"][key]["Query_seq"] = sequences[hit["Query_ORFID"]]

    
    def _format_results(self, results_data, incremental_output):
        """
        Process results into their final CSV format and write them to disk.
        """
        # Remove temporary hit counter tag
        for contig in results_data:
            for neighborhood in results_data[contig]:
                del results_data[contig][neighborhood]["new_hit_count"]
        
        filename = "{}_results.csv".format(self.job_id) if self.job_id is not None else "gene_finder_results.csv"
        if self.output_directory is not None and os.path.exists(self.output_directory):
            filename = os.path.join(self.output_directory, filename)
        csv_writer = CSVWriter(results_data, filename)
        mode = "a" if incremental_output and self._appending_results else "w"
        csv_writer.to_csv(self.data_path, mode)
        self._appending_results = True
    
    
    def _record_all_hits(self, all_hit_data):
        """
        Write all hits from one or more contigs to disk (in json format).
        """
        filename = "{}_hits.json".format(self.job_id) if self.job_id is not None else "gene_finder_hits.json"
        if self.output_directory is not None and os.path.exists(self.output_directory):
            filename = os.path.join(self.output_directory, filename)
        mode = "a" if self._appending_hits else "w"
        # looping allows the json-like data to be written incrementally, which is 
        # admittedly not really what json is intended for
        # we'll do this regardless of whether `incremental_output` is true so
        # that the file format is the same in both cases
        with open(filename, mode) as f:
            for contig in all_hit_data:
                final_candidate_count = self._final_candidate_count(contig)
                contig_data = {contig: {"final_candidate_count": final_candidate_count, "hits": all_hit_data[contig]}}
                json.dump(contig_data, f)
                f.write("\n")
        self._appending_hits = True


    def _final_candidate_count(self, contig_id):
        """
        Returns the number of candidate systems (per contig) that were found.
        """
        return len(self._results[contig_id])
    
    
    def _write_checkpoint_file(self, contig_id):
        """
        Write the ID of the current contig (i.e the contig currently being processed)
        to disk. This gets overwritten with each new contig ID, and removed once the
        job completes successfully. The point is that if the user is running gene 
        finder in "incremental" mode and the job fails unexpectedly, the ID can be
        used to re-start gene finder from where it left off.
        """
        filename = "{}_checkpoint.txt".format(self.job_id) if self.job_id is not None else "gene_finder_checkpoint.txt"
        if self.output_directory is not None and os.path.exists(self.output_directory):
            filename = os.path.join(self.output_directory, filename)
        with open(filename, "w") as f:
            f.write("{},{}".format(contig_id, self.data_path))
    

    def _remove_checkpoint_file(self, incremental_output):
        """
        Remove the contig ID checkpoint file. This is only called after the entire job
        has been successfully completed.
        """
        filename = "{}_checkpoint.txt".format(self.job_id) if self.job_id is not None else "gene_finder_checkpoint.txt"
        if self.output_directory is not None and os.path.exists(self.output_directory):
            filename = os.path.join(self.output_directory, filename)
        # check that the checkpoint file was actually created during this run
        if os.path.exists(filename) and incremental_output:
            os.remove(filename)
    

    def run(self, data, job_id=None, output_directory=None, min_prot_len=60, 
            span=10000, record_all_hits=False, incremental_output=False,
            starting_contig=None, gzip=False) -> dict:
        """
        Execute each step in the pipeline, in the order they were added.

        Args:
            data (str): Path to the input data file. Can be a single-
                or multi-sequence file in fasta format.
            job_id (str, optional): A unique ID to prefix all output
                files. If no ID is given, the string "gene_finder" 
                will be used as the prefix. In any case, results from
                the pipeline are written to the file <prefix>_results.csv.
            output_directory (str, optional): The directory to write
                output data files to. If no directory is given then the current
                (working) directory is used.
            min_prot_len (int, optional): Minimum ORF length (aa).
                Default is 60.
            span (int, optional): Length (nt) upsteam and downstream
                of each seed hit to keep. Defines the aproximate size
                of the genomic neighborhoods that will be used as the
                search space after the seed step.
            record_all_hits (bool, optional): Write data about all genes found 
                (even discarded ones) to the file <job_id>_hits.json,
                grouped by contig. Note that this contains much of the same
                information as is in the results CSV file; nevertheless, it
                may be useful for analysis or troubleshooting a search.
            incremental_output (bool, optional): Write results to disk
                after each contig is processed. Using this option also creates a 
                checkpoint file that gives the ID of the contig that is currently 
                being processed; if the job finishes successfully, this file will 
                be automatically cleaned up. This feature is especially useful 
                for long-running jobs. 
            starting_contig (bool, optional): The sequence identifier of
                the contig where the run should begin. In other words, 
                skip over records in the input file until
                the specified contig is reached, and then run the pipeline
                as normal. This is usually used in conjunction with 
                ``incremental_output``.
            gzip (bool, optional): Was this file compressed with gzip? 

        Returns:   
            dict: Candidate systems, grouped by contig id and genomic location.
        """
        self._reset_results()
        self.data_path = data
        self.min_prot_len = min_prot_len
        self.span = span
        self.job_id = job_id
        self.output_directory = output_directory

        data_handle = self._open_data(self.data_path, gzip)
        for record in SeqIO.parse(data_handle, "fasta"):
            if starting_contig is not None:
                if record.id == starting_contig:
                    starting_contig = None
                else:
                    continue
            # if this run was seeded with a particular contig id, skip all other contigs
            # in the file if they exit
            # saves a bit of time, since there is overhead associated with getting things set up
            if isinstance(self._steps[0], SeedWithCoordinatesStep):
                if self._steps[0].contig_id is not None and self._steps[0].contig_id != record.id:
                    continue
            # clear out any data that may be leftover from processing a previous
            # contig and get all ORFs in this contig
            contig_id, contig_len, contig_path = self._setup_run(record)
            if incremental_output:
                self._write_checkpoint_file(contig_id)
            if self._all_orfs is None:
                # No ORFs were identified in the contig (probably because it's small),
                # so we just continue on to the next contig if it exists
                self._all_hits[contig_id] = {}
                self._results[contig_id] = {}
                if incremental_output and record_all_hits:
                    self._record_all_hits({contig_id: {}})
                continue
            for step in self._steps:
                if isinstance(step, SeedStep):
                    #print("Begin seed step: {}".format(step.name))
                    step.execute(self._all_orfs, self.span, contig_len)
                    self._get_orfs_in_neighborhood(step.neighborhood_ranges, 
                                                    contig_path,
                                                    contig_id)
                    self._results_init(step.neighborhood_ranges)
                    self._results_update(step.hits, min_prot_count=0)
                    neighborhood_orfs = concatenate(self._working_dir.name, 
                                                    self._neighborhood_orfs.values())

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
                    step.execute(neighborhood_orfs)
                    self._results_update(step.hits, min_prot_count=step.min_prot_count)
                    neighborhood_orfs = concatenate(self._working_dir.name, 
                                                    self._neighborhood_orfs.values())
                    #print("Filtering for {} genes complete".format(step.name))
                
                elif isinstance(step, CrisprStep):
                    #print("Begin CRISPR array search")
                    step.execute(contig_path)
                    self._results_update_crispr(step.hits)
                    #print("CRISPR array search complete")
                        
                elif isinstance(step, BlastnStep):
                    step.execute(contig_path)
                    self._results_update_nucl(step.hits)

                else:
                    #print("Begin blast step: {}".format(step.name))
                    step.execute(neighborhood_orfs)
                    self._results_update(step.hits, min_prot_count=0)
                    #print("Blast step complete")
                
                if len(step.hits) != 0:
                    self._all_hits[contig_id][step.search_tool.step_id] = step.hits
                # check if all candidates have been eliminated (if so, no need to
                # continue searching)
                if len(self._neighborhood_orfs) == 0:
                    #print("No putative neighborhoods remain - terminating run")
                    self._results[contig_id] = {}
                    break
                self._update_output_sequences()
                self._results[contig_id] = self._working_results

            # The current contig is finished processing; if outputing results incrementally,
            # we'll dump data about this contig to disk now 
            if incremental_output:
                self._format_results({contig_id: self._results[contig_id]}, True)
                if record_all_hits:
                    self._record_all_hits({contig_id: self._all_hits[contig_id]})

        # all contigs have finished processing; dump all data to disk if we haven't already
        if not incremental_output:
            self._format_results(self._results, False)
            if record_all_hits:
                self._record_all_hits(self._all_hits)
        
        self._remove_checkpoint_file(incremental_output)
        data_handle.close() # make sure to close the input data file object
        return self._results
