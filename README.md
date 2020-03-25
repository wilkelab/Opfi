# crisposon

## Requirements

Crisposon uses NCBI BLAST+ applications to perform homology searches, and a slightly modified version of [piler-cr](https://www.drive5.com/pilercr/) to identify CRISPR arrays. 

You can find NCBI BLAST+ install info [here](https://www.ncbi.nlm.nih.gov/books/NBK279671/).

Instructions for building piler-cr are detailed below.

## Installation

First grab your own fork of crisposon and clone it. Then (from the main project directory) you can do:
```
pip3 install .
pip3 install -r requirements.txt
```
This should install crisposon and its only python dependency, biopython. 

The last thing you need to do is compile the custom piler-cr source code, which is in lib/pilercr1.06.
```
cd lib/pilercr1.06
make
```
Add the pilercr executable to your PATH you're good to go.

## Using crisposon

The crisposon api is centered around the Pipeline class. Typical usage looks something like this:
```python
from crisposon.pipeline import Pipeline

p = Pipeline("my_genome.fasta", name="V_crassostrea", min_prot_len=30, span=10000)
p.add_seed_step(db="data/blast_databases/tnsAB", name="tnsAB", e_val=0.001, type="PROT")
p.add_filter_step(db="data/blast_databases/cas", name="cas", e_val=0.001, type="PROT")
p.add_blast_step(db="data/blast_databases/tnsCD", name="tnsCD", e_val=0.001, type="PROT")
p.add_crispr_step()

results = p.run()
```
The docstrings in pipeline.py provide more details about pipeline usage, paramaters, etc.
Also, I tried to clean up some of the scripts I've been using to run jobs - they can be found under extras/
