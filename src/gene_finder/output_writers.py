import csv

# Header values in old pipeline output
# CSVWriter no longer writes headers, but this list
# is still needed for operon_analyzer to handle
# old output formats
FIELDNAMES = ["Contig", "Locus coordinates", "Feature", 
              "Feature coordinates", "Query ORFID",
              "Strand", "Hit accession", "Hit e-val", 
              "Description", "Query/Consensus sequence"]


class CSVWriter:
    """
    Write pipeline output to a CSV formatted file.
    """
    def __init__(self, results, outfile):
        """
        Args:
            results (dict): All candidates identified by the pipeline.
            outfile (string): The file to write output to.
        """
        self.results = results
        self.id = None
        self.project_id = None
        self.outfile = outfile
    

    def _ret_fieldnames(self):
        # is this still needed?
        return FIELDNAMES
    

    def _get_nucleotide_row(self, neighborhood, hit):
        """
        Format a row representing a hit identified with nucleotide BLAST.
        """
        row = [""] * 22
        row[0] = self.id
        row[1] = "{}..{}".format(neighborhood["Loc_start-pos"],
                                neighborhood["Loc_end-pos"])
        row[2] = hit["Hit_name"]
        row[3] = "{}..{}".format(hit["Query_start-pos"],
                                hit["Query_end-pos"])
        row[4] = ''
        row[5] = hit["Strand"]
        row[6] = hit["Hit_accession"]
        row[7] = hit["Hit_e-val"]
        row[8] = hit["Hit_description"]
        row[9] = hit["Query_seq"]
        row[10] = hit["Bitscore"]
        row[11] = hit["Raw_score"]
        row[12] = hit["Alignment_length"]
        row[13] = hit["Alignment_percent-identical"]
        row[14] = hit["Alignment_num-identical"]
        row[15] = hit["Alignment_mismatches"]
        row[16] = hit["Alignment_num-positive"] if hit["Alignment_num-positive"] is not None else ""
        row[17] = hit["Alignment_num-gapopenings"] 
        row[18] = hit["Alignment_num-gaps"] if hit["Alignment_num-gaps"] is not None else ""
        row[19] = hit["Alignment_percent-pos"] if hit["Alignment_percent-pos"] is not None else ""
        row[20] = hit["Alignment_query-cov"]
        row[21] = self.project_id

        return row
    

    def _get_row(self, neighborhood, hit):
        """
        Format a row representing a hit identified with protein BLAST (or psiBLAST, mmseqs, or 
        Diamond).
        """
        row = [""] * 22
        row[0] = self.id
        row[1] = "{}..{}".format(neighborhood["Loc_start-pos"],
                                neighborhood["Loc_end-pos"])
        row[2] = hit["Hit_name"]
        row[3] = "{}..{}".format(hit["Query_start-pos"],
                                hit["Query_end-pos"])
        row[4] = hit["Query_ORFID"]
        row[5] = hit["Query_ORFID"].split("|")[-1]
        row[6] = hit["Hit_accession"]
        row[7] = hit["Hit_e-val"]
        row[8] = hit["Hit_description"]
        row[9] = hit["Query_seq"]
        row[10] = hit["Bitscore"]
        row[11] = hit["Raw_score"]
        row[12] = hit["Alignment_length"]
        row[13] = hit["Alignment_percent-identical"]
        row[14] = hit["Alignment_num-identical"]
        row[15] = hit["Alignment_mismatches"]
        row[16] = hit["Alignment_num-positive"] if hit["Alignment_num-positive"] is not None else ""
        row[17] = hit["Alignment_num-gapopenings"] 
        row[18] = hit["Alignment_num-gaps"] if hit["Alignment_num-gaps"] is not None else ""
        row[19] = hit["Alignment_percent-pos"] if hit["Alignment_percent-pos"] is not None else ""
        row[20] = hit["Alignment_query-cov"]
        row[21] = self.project_id

        return row
    

    def _format_array_des(self, hit):
        """
        Create a string representation of a CRISPR array hit.
        """
        copy = hit["Copies"]
        repeat = hit["Repeat"]
        spacer = hit["Spacer"]

        array_des = "Copies: {}, Repeat: {}, Spacer: {}".format(copy, repeat, spacer)
        
        return array_des


    def _get_crispr_array_row(self, neighborhood, array):
        """
        Format a row representing a CRISPR array identified using PILER-CR.
        """
        row = [""] * 22
        row[0] = self.id
        row[1] = "{}..{}".format(neighborhood["Loc_start-pos"],
                                                    neighborhood["Loc_end-pos"])
        row[2] = "CRISPR array"
        row[3] = "{}..{}".format(array["Position"],
                                str(int(array["Position"]) + int(array["Length"])))
        row[8] = self._format_array_des(array)
        row[9] = array["Consensus"]
        row[21] = self.project_id

        return row
    

    def _get_rows(self, neighborhood):
        """
        Create a CSV entry for a candidate identified by the pipeline, where each gene/feature 
        is represented by a single row.
        """
        rows = []
        for hit in neighborhood["Hits"]:
            hit_type = neighborhood["Hits"][hit].get('type')
            if hit_type == 'protein':
                rows.append(self._get_row(neighborhood, neighborhood["Hits"][hit]))
            elif hit_type == 'nucleotide':
                rows.append(self._get_nucleotide_row(neighborhood, neighborhood["Hits"][hit]))
            else:
                rows.append(self._get_crispr_array_row(neighborhood, 
                                                        neighborhood["Hits"][hit]))
        
        return rows


    def to_csv(self, project_id, mode):
        """
        Write candidate CSV data to disk.

        Args:
            project_id (string): The filename/path for the input data used by this run 
                of the pipeline. 
            mode (string): The mode for opening the output file for writing. Usually this
                will be "w" (for create/write mode), but if the pipeline is being run 
                incrementally then results may need to be appended to the output file 
                periodically.
        """
        with open(self.outfile, mode, newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            self.project_id = project_id
            for contig in self.results:
                self.id = contig
                for neighborhood in self.results[contig]:
                    rows = self._get_rows(self.results[contig][neighborhood])
                    for row in rows:
                        writer.writerow(row)
