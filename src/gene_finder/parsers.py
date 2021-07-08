import os, json, csv


def parse_blastn_output(tsv, step_id, parse_descriptions=True):
    """
    Parse output from a search step (in BLAST tabular format).

    Expects a tsv file with the following format fields: 
        qseqid sseqid stitle evalue \
        bitscore score length pident \
        nident mismatch positive gapopen \
        gaps ppos qcovhsp qseq qstart qend sstrand

    Returns a dictionary of best hits for each query that had a 
    hit, where "best" means the lowest e-value score.
    """
    try:
        hits = {}
        with open(tsv, newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")

            for row in reader:
                if _keep_row(row, hits):
                    hit_dic = {}
                    hit_dic['type'] = 'nucleotide'

                    # information about the reference protein
                    hit_acc = row[1].strip()
                    hit_dic["Hit_accession"] = hit_acc if len(hit_acc) > 0 else "No accession/ID found for reference"
                    hit_def = row[2].split()
                    if parse_descriptions and len(hit_def) >= 2:
                        hit_dic["Hit_name"] = hit_def[1]
                        hit_dic["Hit_description"] = " ".join(hit_def)
                    else:
                        hit_def = " ".join(hit_def)
                        hit_dic["Hit_name"] = hit_def if len(hit_def) > 0 else "No gene name found for reference"
                        hit_dic["Hit_description"] = hit_def if len(hit_def) > 0 else "No gene description found for reference"
                    
                    # alignment statistics
                    hit_dic["Hit_e-val"] = row[3]
                    hit_dic["Bitscore"] = row[4]
                    hit_dic["Raw_score"] = row[5]
                    hit_dic["Alignment_length"] = row[6]
                    hit_dic["Alignment_percent-identical"] = row[7]
                    hit_dic["Alignment_num-identical"] = row[8]
                    hit_dic["Alignment_mismatches"] = row[9]
                    hit_dic["Alignment_num-positive"] = row[10]
                    hit_dic["Alignment_num-gapopenings"] = row[11]
                    hit_dic["Alignment_num-gaps"] = row[12]
                    hit_dic["Alignment_percent-pos"] = row[13]
                    hit_dic["Alignment_query-cov"] = row[14]
                    hit_dic["Query_seq"] = row[15]
                    hit_dic["Query_start-pos"] = row[16]
                    hit_dic["Query_end-pos"] = row[17]
                    hit_dic["Strand"] = 1 if row[18] == 'plus' else -1
                    hit_name = row[0]
                    hits[hit_name] = hit_dic

        hits = _reformat_hit_ids(hits, step_id)
        return hits
    
    except csv.Error:
        return {}


def parse_search_output(tsv, step_id, search_type, parse_descriptions=True):
    """
    Parse output from a search step (in BLAST tabular format).

    Expects a tsv file with the following format
    fields: qseqid sseqid stitle evalue \
        bitscore score length pident \
        nident mismatch positive (blast/diamond only) \
        gapopen gaps (blast/diamond only) \
        ppos (blast/diamond only) qcovs

    Returns a dictionary of best hits for each query that had a 
    hit, where "best" means the lowest e-value score.
    """
    try:
        hits = {}
        with open(tsv, newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")

            for row in reader:
                if _keep_row(row, hits):
                    hit_dic = {}
                    local_query_id = row[0].split()[0]
                    hit_dic['type'] = 'protein'
                    hit_dic["Query_ORFID"] = local_query_id

                    # get query start/stop pos (nt) from ORF id
                    local_query_id = local_query_id.split("|")
                    hit_dic["Query_start-pos"] = local_query_id[1]
                    hit_dic["Query_end-pos"] = local_query_id[2]

                    # information about the reference protein
                    hit_acc = row[1].strip()
                    hit_dic["Hit_accession"] = hit_acc if len(hit_acc) > 0 else "No accession/ID found for reference"
                    hit_def = row[2].split()
                    if parse_descriptions and len(hit_def) >= 2:
                        hit_dic["Hit_name"] = hit_def[1]
                        hit_dic["Hit_description"] = " ".join(hit_def)
                    else:
                        hit_def = " ".join(hit_def)
                        hit_dic["Hit_name"] = hit_def if len(hit_def) > 0 else "No gene name found for reference"
                        hit_dic["Hit_description"] = hit_def if len(hit_def) > 0 else "No gene description found for reference"
                    
                    # alignment statistics
                    hit_dic["Hit_e-val"] = row[3]
                    hit_dic["Bitscore"] = row[4]
                    hit_dic["Raw_score"] = row[5]
                    hit_dic["Alignment_length"] = row[6]
                    hit_dic["Alignment_percent-identical"] = row[7]
                    hit_dic["Alignment_num-identical"] = row[8]
                    hit_dic["Alignment_mismatches"] = row[9]

                    if search_type in ("blast", "diamond"):
                        hit_dic["Alignment_num-positive"] = row[10]
                        hit_dic["Alignment_num-gapopenings"] = row[11]
                        hit_dic["Alignment_num-gaps"] = row[12]
                        hit_dic["Alignment_percent-pos"] = row[13]
                        hit_dic["Alignment_query-cov"] = row[14]
                        hit_dic["Query_seq"] = row[15]
                    
                    else: # mmseqs doesn't output some of these
                        hit_dic["Alignment_num-positive"] = None
                        hit_dic["Alignment_num-gapopenings"] = row[10]
                        hit_dic["Alignment_num-gaps"] = None
                        hit_dic["Alignment_percent-pos"] = None
                        hit_dic["Alignment_query-cov"] = row[11]
                        hit_dic["Query_seq"] = row[12]

                    hit_name = row[0]
                    hits[hit_name] = hit_dic

        hits = _reformat_hit_ids(hits, step_id)
        return hits
    
    except csv.Error:
        return {}


def _is_int(string):
    """
    Check if a string can be cast to an int.
    """
    try: 
        int(string)
        return True
    except ValueError:
        return False


def _get_column_lables(line):
    """
    Extract the following column lables from 
    pilercr "SUMMARY BY POSITION" table: 
    "Position", "Length", "Copies", "Repeat", "Spacer", "Strand",
    "Consensus".
    """
    lables = line.split()
    lables_to_remove = ["Array", "Sequence", "#", "+"]
    lables = [lable for lable in lables if lable not in lables_to_remove]
    lables.insert(5, "Strand")
    
    return lables
 

def _get_array_info(line, array_num):
    """Extract the following information for an array
    in pilercr "SUMMARY BY SIMILARITY" table: 
    "Position", "Length", "Copies", "Repeat", "Spacer", "Strand",
    "Consensus".
    """
    line = line[25:]
    array_info = line.split()
    array_id = "Array_{}".format(array_num)
    
    return array_id, array_info


def parse_pilercr_summary(pilercr_out):
    """
    Parse pilercr output in the "summary by position" section only.

    Returns a dictionary of dictionaries, containing
    key-value pairs of infomation associated with
    each array identified by piler.

    Example:

        >>> arrays = parse_pilercr_summary("vcrass_pilercr_out")
        >>> print(json.dumps(arrays, indent=4))
        {
            "array_1": {
                "Position": "30001",
                "Length": "148",
                "Copies": "3",
                "Repeat": "28",
                "Spacer": "32",
                "Consensus": "GTGAACTGCCGAATAGGTAGCTGATAAT"
            }
        }

    """    
    headers_reached = False
    array_info_start = False
    lables = []
    arrays = {}
    array_num = 0

    with open(pilercr_out) as f:
        for line in f:

            if headers_reached and not line.isspace():
                if line.strip().split()[0] == "Array":
                    lables = _get_column_lables(line)
                    array_info_start = True
                    headers_reached = False

            elif array_info_start and not line.isspace():
                if _is_int(line.strip().split()[0]):
                    array_id, array_info = _get_array_info(line, array_num)
                    arrays[array_id] = {}
                    array_num += 1
                    for lable, value in zip(lables, array_info):
                        arrays[array_id][lable] = value

            if line.strip() == "SUMMARY BY SIMILARITY":
                headers_reached = True

            if line.strip() == "SUMMARY BY POSITION":
                array_info_start = False

    return arrays


def _keep_row(row, hits):
    """
    Determine if given row represents the best hit for this
    query so far by comparing e-values. 
    """
    if row[0] in hits:
        if float(row[3]) < float(hits[row[0]]["Hit_e-val"]):
            return True
        else:
            return False
    else:
        return True


def _reformat_hit_ids(hits, step_id):
    """
    For consistency, copy over hit info from a dict where keys are the
    utility-specific query IDs to a dict where keys are in the 
    format created by the other gene finder parsers.
    """
    re_hits = {}
    hit_num = 0
    for hit in hits.values():
        hit_id = "{}{}_hit-{}".format(step_id[:1].upper(), step_id[1:], str(hit_num))
        re_hits[hit_id] = hit
        hit_num += 1
    
    return re_hits
