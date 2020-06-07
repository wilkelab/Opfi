import csv

FIELDNAMES = ["Contig", "Locus coordinates", "Feature", 
              "Feature coordinates", "Query ORFID",
              "Strand", "Hit accession", "Hit e-val", 
              "Description", "Query/Consensus sequence"]


class CSVWriter:

    def __init__(self, results, outfile):

        self.results = results
        self.id = None
        self.project_id = None
        self.outfile = outfile
    
    def _ret_fieldnames(self):

        return FIELDNAMES
    
    def _get_row(self, neighborhood, hit):
        
        row = [""] * 21
        row[0] = self.project_id
        row[1] = self.id
        row[2] = "{}..{}".format(neighborhood["Loc_start-pos"],
                                neighborhood["Loc_end-pos"])
        row[3] = hit["Hit_name"]
        row[4] = hit["Query_ORFID"]
        row[5] = "{}..{}".format(hit["Query_start-pos"],
                                hit["Query_end-pos"])
        row[6] = hit["Query_ORFID"].split("|")[-1]
        row[7] = hit["Hit_accession"]
        row[8] = hit["Hit_description"]
        row[9] = hit["Hit_e-val"]
        row[10] = hit["Bitscore"]
        row[11] = hit["Raw_score"]
        row[12] = hit["Alignment_length"]
        row[13] = hit["Alignment_percent-identical"]
        row[14] = hit["Alignment_num-identical"]
        row[15] = hit["Alignment_mismatches"]
        row[16] = hit["Alignment_num-positive"]
        row[17] = hit["Alignment_num-gapopenings"]
        row[18] = hit["Alignment_num-gaps"]
        row[19] = hit["Alignment_percent-pos"]
        row[20] = hit["Alignement_query-cov"]

        return row
    
    def _format_array_des(self, hit):

        copy = hit["Copies"]
        repeat = hit["Repeat"]
        spacer = hit["Spacer"]

        array_des = "Copies: {}, Repeat: {}, Spacer: {}".format(copy, repeat, spacer)
        
        return array_des

    def _get_crispr_array_row(self, neighborhood, array):
        
        row = [""] * 21
        row[0] = self.project_id
        row[1] = self.id
        row[2] = "{}..{}".format(neighborhood["Loc_start-pos"],
                                                    neighborhood["Loc_end-pos"])
        row[3] = "CRISPR array"
        row[5] = "{}..{}".format(array["Position"],
                                str(int(array["Position"]) + int(array["Length"])))
        row[8] = self._format_array_des(array)

        return row
    
    def _get_rows(self, neighborhood):

        rows = []
        for hit in neighborhood["Hits"]:
            if "Query_ORFID" in neighborhood["Hits"][hit]:
                rows.append(self._get_row(neighborhood, neighborhood["Hits"][hit]))
            else:
                rows.append(self._get_crispr_array_row(neighborhood, 
                                                        neighborhood["Hits"][hit]))
        
        return rows

    def to_csv(self, project_id):

        with open(self.outfile, 'w', newline='') as csvfile:
            fieldnames = self._ret_fieldnames()
            writer = csv.writer(csvfile, fieldnames=fieldnames)
            
            self.project_id = project_id
            for contig in self.results:
                self.id = contig
                for neighborhood in self.results[contig]:
                    rows = self._get_rows(self.results[contig][neighborhood])
                    for row in rows:
                        writer.writerow(row)
