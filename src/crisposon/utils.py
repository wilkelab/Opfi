import os, shutil

def concatenate(in_dir, file_names):
    out = os.path.join(in_dir, "merged_input.fasta")
    if os.path.exists(out):
        os.remove(out)
    
    with open(out,"wb") as outfile:
        for name in file_names:
            with open(name,"rb") as infile:
                shutil.copyfileobj(infile, outfile)
    
    return out

def get_neighborhood_ranges(hits, span=20000):

    hit_coords = []
    keys_sorted = sorted(hits, key=lambda k: min(int(hits[k]["q_start"]), int(hits[k]["q_stop"])))
    
    for key in keys_sorted:
        start = int(hits[key]["q_start"])
        stop = int(hits[key]["q_stop"])
        lower = min(start, stop)
        upper = max(start, stop)
        hit_coords.append((lower, upper))
    
    ranges = []
    first_coord = hit_coords.pop(0)
    lower = max((first_coord[0] - span), 0)
    upper = first_coord[1] + span
    ranges.append((lower, upper))
    
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

    