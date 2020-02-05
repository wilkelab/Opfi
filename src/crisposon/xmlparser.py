import os
import xml.etree.ElementTree as ET

def parse_blast(blast_xml, blast_id):
    
    tree = ET.parse(blast_xml)
    blast_output = tree.getroot()
    
    hits = {}
    hit_num = 0
    for iteration in blast_output.iter('Iteration'):
        hit = iteration.find('Iteration_hits').find('Hit')
        
        if hit is not None:
            hit_dic = {}
            
            query = iteration.find("Iteration_query-def").text
            query = query.split().pop(0)
            hit_dic["q_id"] = query

            query = query.split("|")
            hit_dic["q_start"] = query[1]
            hit_dic["q_stop"] = query[2]

            hit_def = hit.find("Hit_def").text.split()
            hit_dic["common_name"] = hit_def.pop(1)
            hit_dic["ref_acc"] = hit_def.pop(0)
            hit_dic["ref_des"] = " ".join(hit_def)

            e_val = hit.find("Hit_hsps").find("Hsp").find("Hsp_evalue").text
            hit_dic["e_val"] = e_val

            hit_len = hit.find("Hit_len").text
            hit_dic["hit_len"] = hit_len

            hit_name = blast_id + "_hit_" + str(hit_num)
            hits[hit_name] = hit_dic
            hit_num += 1

    return hits
