import os, json, csv
import xml.etree.ElementTree as ET

def parse_blast_xml(blast_xml, blast_id):
    """Parse blast xml output generated from 
    running a blast step in the pipeline. 

    Queries (ORFs) for which no alignments were returned
    are excluded from the output.

    Returns a dictionary of dictionaries if hits
    were found; otherwise, returns None.

    Args:
        blast_xml(str): Path to blast output.
        blast_id(str, optional): A unique ID for the
            blast step.
    """
    try:
        tree = ET.parse(blast_xml)
        blast_output = tree.getroot()
        
        hits = {}
        hit_num = 0
        for iteration in blast_output.iter('Iteration'):
            hit = iteration.find('Iteration_hits').find('Hit')
            
            if hit is not None:
                hit_dic = {}
                
                # get the ORF id
                query = iteration.find("Iteration_query-def").text
                query = query.split().pop(0)
                hit_dic["Query_ORFID"] = query

                # get query start/stop pos (nt) from ORF id
                query = query.split("|")
                hit_dic["Query_start-pos"] = query[1]
                hit_dic["Query_end-pos"] = query[2]

                # information about the reference protein
                hit_def = hit.find("Hit_def").text.split()
                hit_dic["Hit_name"] = hit_def.pop(1)
                hit_dic["Hit_accession"] = hit_def.pop(0)
                hit_dic["Hit_description"] = " ".join(hit_def)

                # e val for this hsp
                e_val = hit.find("Hit_hsps").find("Hsp").find("Hsp_evalue").text
                hit_dic["Hit_e-val"] = e_val

                # Sequence of the translated ORF used as the query
                query_sequence = hit.find("Hit_hsps").find("Hsp").find("Hsp_qseq").text
                hit_dic["Query_sequence"] = query_sequence

                # Capitalize first letter only of blast_id (set by pipeline.add_step() name param)
                hit_name = "{}{}_hit-{}".format(blast_id[:1].upper(), blast_id[1:], str(hit_num))
                hits[hit_name] = hit_dic
                hit_num += 1

        return hits
    
    except ET.ParseError:

        return {}

def _is_int(string):
    """Check if a string can be cast to an int.

    Helps parse_pilercr_summary determine if 
    a line contains array information.
    """
    try: 
        int(string)
        return True
    except ValueError:
        return False

def _get_column_lables(line):
    """Extract the following column lables from 
    pilercr "SUMMARY BY POSITION" table: "Position",
    "Length", "Copies", "Repeat", "Spacer", "Strand",
     "Consensus".

    Used by parse_pilercr_summary.
    """
    lables = line.split()
    lables_to_remove = ["Array", "Sequence", "#", "+"]
    lables = [lable for lable in lables if lable not in lables_to_remove]
    lables.insert(5, "Strand")
    
    return lables
 
def _get_array_info(line, array_num):
    """Extract the following information for an array
    in pilercr "SUMMARY BY SIMILARITY" table: "Position",
    "Length", "Copies", "Repeat", "Spacer", "Strand",
    "Consensus".

    Used by parse_pilercr_summary.
    """
    line = line[25:]
    array_info = line.split()
    array_id = "Array_{}".format(array_num)
    
    return array_id, array_info

def parse_pilercr_summary(pilercr_out):
    """Parse pilercr output.
    
    "summary by position" section only.

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

def _keep_row(row, hits, queryid):
    """Parser helper function.

    Determine if given row represents the best hit for this
    query so far. 
    """

    if row[queryid] in hits:
        if float(row["evalue"]) < float(hits[row[queryid]]["Hit_e-val"]):
            return True
        else:
            return False
    else:
        return True

def _reformat_hit_ids(hits, step_id):
    """Parser helper function.

    For consistency, copy over hit info from a dict where keys are the
    utility-specific query IDs to a dict where keys are in the 
    format created by the other crisposon parsers.
    """

    re_hits = {}
    hit_num = 0
    for hit in hits.values():
        hit_id = "{}{}_hit-{}".format(step_id[:1].upper(), step_id[1:], str(hit_num))
        re_hits[hit_id] = hit
        hit_num += 1
    
    return re_hits

def parse_mmseqs(mmseqs_tsv, step_id, fields):
    """Parse mmseqs output (in blast tabular format).

    Expects an mmseqs output tsv file with the following 
    fields: "query,target,evalue,qseq,qheader,theader,qcov,tset"

    Returns a dictionary of best hits for each query that had a 
    hit, where "best" means the lowest e-value score.
    """

    try:
        hits = {}
        with open(mmseqs_tsv, newline='') as tsvfile:
            reader = csv.DictReader(tsvfile, fieldnames=fields, delimiter="\t")

            for row in reader:
                if _keep_row(row, hits, "query"):
                    hit_dic = {}
                    local_query_id = row["qheader"].split()[0]
                    hit_dic["Query_ORFID"] = local_query_id

                    # get query start/stop pos (nt) from ORF id
                    local_query_id = local_query_id.split("|")
                    hit_dic["Query_start-pos"] = local_query_id[1]
                    hit_dic["Query_end-pos"] = local_query_id[2]

                    # information about the reference protein
                    hit_dic["Hit_description"] = row["theader"]
                    hit_dic["Hit_name"] = step_id
                    hit_dic["Hit_accession"] = row["target"]
                    
                    # e val for this hsp
                    hit_dic["Hit_e-val"] = row["evalue"]

                    # Sequence of the translated ORF used as the query
                    hit_dic["Query_sequence"] = row["qseq"]

                    hit_name = row["query"]
                    hits[hit_name] = hit_dic

        hits = _reformat_hit_ids(hits, step_id)
        return hits
    
    except csv.Error:
        return {}

def parse_diamond(diamond_tsv, step_id, fields):
    """Parse diamond output (in blast tabular format).

    Expects a diamond output tsv file with the following 
    fields: qseqid sseqid full_qseq evalue stitle

    Returns a dictionary of best hits for each query that had a 
    hit, where "best" means the lowest e-value score.
    """

    try:
        hits = {}
        with open(diamond_tsv, newline='') as tsvfile:
            reader = csv.DictReader(tsvfile, fieldnames=fields, delimiter="\t")

            for row in reader:
                if _keep_row(row, hits, "qseqid"):
                    hit_dic = {}
                    local_query_id = row["qseqid"]
                    hit_dic["Query_ORFID"] = local_query_id

                    # get query start/stop pos (nt) from ORF id
                    local_query_id = local_query_id.split("|")
                    hit_dic["Query_start-pos"] = local_query_id[1]
                    hit_dic["Query_end-pos"] = local_query_id[2]

                    # information about the reference protein
                    hit_dic["Hit_description"] = row["stitle"]
                    hit_dic["Hit_name"] = step_id
                    hit_dic["Hit_accession"] = row["sseqid"]
                    
                    # e val for this hsp
                    hit_dic["Hit_e-val"] = row["evalue"]

                    # Sequence of the translated ORF used as the query
                    hit_dic["Query_sequence"] = row["full_qseq"]

                    hit_name = row["qseqid"]
                    hits[hit_name] = hit_dic

        hits = _reformat_hit_ids(hits, step_id)
        return hits
    
    except csv.Error:
        return {}