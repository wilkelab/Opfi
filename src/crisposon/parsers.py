import os, json
import xml.etree.ElementTree as ET

def parse_blast_xml(blast_xml, blast_id):
    """Parse blast xml output generated from 
    running a blast step in the pipeline. 

    Queries (ORFs) for which no alignments were returned
    are excluded from the output.

    Returns a dictionary of dictionaries if hits 
    were found; otherwise, returns an empty dictionary.

    Args:
        blast_xml(str): Path to blast output.
        blast_id(str, optional): A unique ID for the
            blast step.
    """
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
            hit_dic["ORFID"] = query

            # get query start/stop pos (nt) from ORF id
            query = query.split("|")
            hit_dic["Start"] = query[1]
            hit_dic["Stop"] = query[2]

            # information about the reference protein
            hit_def = hit.find("Hit_def").text.split()
            hit_dic["Name"] = hit_def.pop(1)
            hit_dic["Accession"] = hit_def.pop(0)
            hit_dic["Description"] = " ".join(hit_def)

            # e val for this hsp
            e_val = hit.find("Hit_hsps").find("Hsp").find("Hsp_evalue").text
            hit_dic["Expect"] = e_val

            # Sequence of the translated ORF used as the query
            query_sequence = hit.find("Hit_hsps").find("Hsp").find("Hsp_qseq").text
            hit_dic["QuerySeq"] = query_sequence

            hit_name = blast_id + "_hit_" + str(hit_num)
            hits[hit_name] = hit_dic
            hit_num += 1

    return hits

def _get_field_lables(line):
    """Extract the following field lables from 
    pilercr "SUMMARY BY POSITION" table: "Position",
    "Length", "Copies", "Repeat", "Spacer", 
    "Consensus".

    Used by parse_pilercr_summary.
    """
    fields = line.split()
    fields_to_remove = ["Array", "Sequence", "#", "Distance"]
    fields = [field for field in fields if field not in fields_to_remove]
    
    return fields

def _get_array_info(line):
    """Extract the following information for an array
    in pilercr "SUMMARY BY POSITION" table: "Position",
    "Length", "Copies", "Repeat", "Spacer", 
    "Consensus".

    Used by parse_pilercr_summary.
    """
    array_info = line.split("  ")
    array_info = map(str.strip, array_info)
    array_info = [info for info in array_info if info != ""]
    if len(array_info) == 9:
        array_info.pop(7)
    array_info.pop(1)
    array_num = array_info.pop(0)
    array_id = "array_{}".format(array_num)
    
    return array_id, array_info

def parse_pilercr_summary(pilercr_out):
    """Parse pilercr output - "summary
    by position" section only.

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
    summary_reached = False
    array_info_start = False
    fields = []
    arrays = {}

    with open(pilercr_out) as f:
        for line in f:
            line = line.strip()

            if summary_reached and not len(line.split()) == 0:
                if line.split()[0] == "Array":
                    fields = _get_field_lables(line)
                    array_info_start = True
                    summary_reached = False

            elif array_info_start and not len(line.split()) == 0:
                if line.split()[0] != "==========":
                    array_id, array_info = _get_array_info(line)
                    arrays[array_id] = {}
                    for field, value in zip(fields, array_info):
                        arrays[array_id][field] = value

            elif line == "SUMMARY BY POSITION":
                summary_reached = True

    return arrays

if __name__ == "__main__":
    #pilercr_out = "/home/alexis/Projects/CRISPR-Transposons/out/pilercr/vcrass"
    #pilercr_out = "/home/alexis/Projects/CRISPR-Transposons/out/pilercr/a_brierleyi"

    #arrays = parse_pilercr_summary(pilercr_out)
    #print(json.dumps(arrays, indent=4))

    blast_out = "/home/alexis/Projects/CRISPR-Transposons/data/tmp/v_crass.xml"
    hits = parse_blast_xml(blast_out, blast_id="v_crass")
    print(json.dumps(hits, indent=4))