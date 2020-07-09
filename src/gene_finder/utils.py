import os, shutil

def concatenate(in_dir, file_names):
    """Concatenate two or more fasta files.

    Merged file is written to the same
    directory as the input files.

    Args:
        in_dir (str): Path to input file directory.
        file_names (list): List of file names.
    """
    out = os.path.join(in_dir, "merged_input.fasta")
    if os.path.exists(out):
        os.remove(out)
    
    with open(out,"wb") as outfile:
        for name in file_names:
            with open(name,"rb") as infile:
                shutil.copyfileobj(infile, outfile)
    
    return out

    
def get_neighborhood_ranges(hits, contig_len, span=20000):
    """Determine the start and end positions of genomic
    neighborhoods surrounding one or more blast hit.

    Neighborhood size is set by the span parameter; specifically,
    this is the number of nucleotides directly to the left and right
    of the hit.

    The position of the hit is given by the corresponding ORF used 
    as the query.

    When two or more hits have overlapping neighborhood regions, 
    the result is merged - i.e, a single neighborhood range is
    returned for these hits.

    Args:
        hits (dict): Parsed blast output. 
        span (int, optional): Number of nucleotides directly to the 
            left and right of the hit to retain. Default is 20000.
    """
    # sort the hit dictionary keys by hit start position
    keys_sorted = sorted(hits, key=lambda k: min(int(hits[k]["Query_start-pos"]), int(hits[k]["Query_end-pos"])))
    
    tmp_ranges = []
    # use sorted keys to construct a temporary list of ranges
    for key in keys_sorted:
        start = int(hits[key]["Query_start-pos"])
        stop = int(hits[key]["Query_end-pos"])
        lower = max(min((start - span), (stop - span)), 0)
        upper = min(max((start + span), (stop + span)), contig_len)
        tmp_ranges.append((lower , upper))
    
    final = []
    # merge overlapping ranges to create final list
    for this_coord in tmp_ranges:
        if len(final) == 0:
            final.append(this_coord)
            continue
        
        last_coord = final[-1]
        if this_coord[0] > last_coord[1]:
            final.append(this_coord)
        elif this_coord[0] <= last_coord[1] and this_coord[1] > last_coord[1]:
            final[-1] = (last_coord[0], this_coord[1])
        else:
            continue

    return final
