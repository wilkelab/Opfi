import os, subprocess, shutil

def concatenate(in_dir, file_names):
    """Merge two or more fasta files."""
    
    path_to_mergefile = os.path.join(in_dir, "merged_input.fasta")

    with open(path_to_mergefile,"wb") as mergefile:
        for name in file_names:
            path_to_infile = os.path.join(in_dir, name)
            with open(path_to_infile,"rb") as infile:
                shutil.copyfileobj(infile, mergefile)
    
    return path_to_mergefile

def reader(path):
    """Determine if the input should be merged."""
    
    merged = False
    if os.path.isfile(path):
        pass
    
    elif os.path.isdir(path):
        files = os.listdir(path)
        if len(files) == 0:
            raise Exception("No files present in input directory")
        
        elif len(files) > 1:
            m_file = concatenate(path, files)
            path = os.path.join(path, m_file)
            merged = True
        
        else:
            path = os.path.join(path, files[0])
    
    else:
        raise Exception("Invalid input")

    return path, merged

def build_blast_db(input, db_dir=None, db_name="blast_db"):
    in_file, merged = reader(input)
    
    if db_dir == None:
        db_dir = os.path.join(input, "blast_db")
    
    if not os.path.exists(db_dir):
        os.mkdir(db_dir)
    
    db = os.path.join(db_dir, db_name)
    cmd = ['makeblastdb', '-in', in_file, '-dbtype', 'prot', '-out', db, '-hash_index']
    subprocess.run(cmd, check=True)

    if merged:
        os.remove(in_file)

if __name__ == "__main__":
    
    ref_dir = "/home/alexis/Projects/CRISPR-Transposons/data/protein_references/cas_uniref"
    db_home_dir = "/home/alexis/Projects/CRISPR-Transposons/data/blast_databases/cas_ALL_sep"

    for fil in os.listdir(ref_dir):

        gene = fil.split("-")[1]
        db_dir = os.path.join(db_home_dir, gene)
        os.mkdir(db_dir)
        build_blast_db(os.path.join(ref_dir, fil), db_dir)
