from Bio import SeqIO
import os, subprocess, shutil

"""def concatenate(in_dir, file_names):
    outstr = "merged_input.fasta"
    
    with open(str(in_dir + "/" + outstr),"wb") as outfile:
        for name in file_names:
            with open(str(in_dir + "/" + name),"rb") as infile:
                shutil.copyfileobj(infile, outfile)
    
    return outstr"""

def concatenate(in_dir, file_names):
    out = os.path.join(in_dir, "merged_input.fasta")
    
    with open(out,"wb") as outfile:
        for name in file_names:
            with open(name,"rb") as infile:
                shutil.copyfileobj(infile, outfile)
    
    return out

def reader(path):
    files = os.listdir(path)
    merged = False

    if len(files) == 0:
        raise Exception("No files present in input directory")
    else:
        if len(files) > 1:
            m_file = concatenate(path, files)
            path = os.path.join(path, m_file)
            merged = True
        else:
            path = os.path.join(path, files[0])
    
    return path, merged

def build_blastp_db(input, db_dir=None, db_name="blast_db"):
    in_file, merged = reader(input)
    
    if db_dir == None:
        db_dir = os.path.join(input, "blast_db")
    
    if not os.path.exists(db_dir):
        os.mkdir(db_dir)
    
    db = os.path.join(db_dir, db_name)
    blastdb_params = ['makeblastdb', '-in', in_file, '-dbtype', 'prot', '-out', db, '-hash_index']
    blastdb_cmd = " ".join(blastdb_params)
    subprocess.call(blastdb_cmd, shell=True)

    if merged:
        os.remove(in_file)


def get_neighborhood_ranges(hits, span=20000):

    hit_coords = []
    keys_sorted = sorted(hits, key=lambda k: min(int(hits[k]["query_start"]), int(hits[k]["query_stop"])))
    
    for key in keys_sorted:
        start = int(hits[key]["query_start"])
        stop = int(hits[key]["query_stop"])
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

    