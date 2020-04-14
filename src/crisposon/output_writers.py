import csv

class CSVWriter:

    def __init__(self, results, contig_id, outfile):

        self.results = results
        self.id = contig_id
        self.outfile = outfile
    
    def _ret_fieldnames(self):

        fieldnames = ["Contig", "Locus coordinates", "Feature", 
                        "Feature coordinates", "Query ORFID",
                        "Strand", "Hit accession", "Hit e-val", 
                        "Description", "Query/Consensus sequence"]
        
        return fieldnames
    
    def _get_row(self, neighborhood, hit):
        
        row = {}
        row["Contig"] = self.id
        row["Locus coordinates"] = "{}..{}".format(neighborhood["Loc_start-pos"],
                                                    neighborhood["Loc_end-pos"])
        row["Feature"] = hit["Hit_name"]
        row["Query ORFID"] = hit["Query_ORFID"]
        row["Feature coordinates"] = "{}..{}".format(hit["Query_start-pos"],
                                                        hit["Query_end-pos"])
        row["Strand"] = hit["Query_ORFID"].split("|")[-1]
        row["Hit accession"] = hit["Hit_accession"]
        row["Hit e-val"] = hit["Hit_e-val"]
        row["Description"] = hit["Hit_description"]
        row["Query/Consensus sequence"] = hit["Query_sequence"]

        return row
    
    def _format_array_des(self, hit):

        copy = hit["Copies"]
        repeat = hit["Repeat"]
        spacer = hit["Spacer"]

        array_des = "Copies: {}, Repeat: {}, Spacer: {}".format(copy, repeat, spacer)
        
        return array_des

    def _get_crispr_array_row(self, neighborhood, array):
        
        row = {}
        row["Contig"] = self.id
        row["Locus coordinates"] = "{}..{}".format(neighborhood["Loc_start-pos"],
                                                    neighborhood["Loc_end-pos"])
        row["Feature"] = "CRISPR array"
        row["Feature coordinates"] = "{}..{}".format(array["Position"],
                                                        str(int(array["Position"]) + int(array["Length"])))
        row["Description"] = self._format_array_des(array)
        row["Query/Consensus sequence"] = array["Consensus"]

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

    def to_csv(self):

        with open(self.outfile, 'w', newline='') as csvfile:
            fieldnames = self._ret_fieldnames()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for neighborhood in self.results:
                rows = self._get_rows(self.results[neighborhood])
                for row in rows:
                    writer.writerow(row)
