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

def get_neighborhood_ranges(hits, span=20000):
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
    hit_coords = []
    # sort the hit dictionary dictionary keys by hit start position
    keys_sorted = sorted(hits, key=lambda k: min(int(hits[k]["q_start"]), int(hits[k]["q_stop"])))
    
    # use sorted keys to construct a sorted list of hit coordinates
    for key in keys_sorted:
        start = int(hits[key]["q_start"])
        stop = int(hits[key]["q_stop"])
        lower = min(start, stop)
        upper = max(start, stop)
        hit_coords.append((lower, upper))
    
    # Construct the first neighborhood 
    ranges = []
    first_coord = hit_coords.pop(0)
    lower = max((first_coord[0] - span), 0)
    upper = first_coord[1] + span
    ranges.append((lower, upper))
    
    # iterate through the remaining hit coordinates, checking if they
    # should be placed in the previous neighborhood or used to
    # construct a new neighborhood
    for coord in hit_coords:
        found_new = False
        last_range = ranges[-1]
        if coord[0] > last_range[1] + span:
            found_new = True
        elif (coord[0] <= last_range[1] + span) and not (coord[1] <= last_range[1]):
            adjusted = (last_range[0], coord[1] + span)
            ranges[-1] = adjusted

        if found_new:
            lower = coord[0] - span
            upper = coord[1] + span
            ranges.append((lower, upper))

    return ranges
    