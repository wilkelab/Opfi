import os, subprocess, shutil

def concatenate(in_dir, file_names):
    outstr = "merged_input.fasta"
    
    with open(str(in_dir + "/" + outstr),"wb") as outfile:
        for name in file_names:
            with open(str(in_dir + "/" + name),"rb") as infile:
                shutil.copyfileobj(infile, outfile)
    
    return outstr

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