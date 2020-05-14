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

p = Pipeline()
p.add_seed_step(db="data/blast_databases/tnsAB", name="tnsAB", e_val=0.001, type="PROT")
p.add_filter_step(db="data/blast_databases/cas", name="cas", e_val=0.001, type="PROT", min_prot_count=2)
p.add_blast_step(db="data/blast_databases/tnsCD", name="tnsCD", e_val=0.001, type="PROT")
p.add_crispr_step()

results = p.run(data="my_genome.fasta", min_prot_len=30, span=10000, outfrmt="CSV", outfile="mygenome.csv")
```
The docstrings in pipeline.py provide more details about pipeline usage, parameters, etc.
Also, I tried to clean up some of the scripts I've been using to run jobs - they can be found under extras/

## Update 05/14

This branch now has some major-ish changes that will eventually be merged into master:

- MMseqs2 and Diamond can now be used as the sequence search tool for the seed, blast, and filter steps
- The pipeline can accept a mult-fasta file as input (no extra work by the user is requried)
- **The pipeline initialization arguments were moved to `pipeline.run()`.** The defaults for `min_prot_len` and `span` are now 60 and 15000 (they were 30 and 10000) respectively.
- **Also, I renamed the `genome` arg to `data` to be more generic (see example above).**
- `name` is no longer an argument of the constructor, nor was it moved to run. Contig id's are now extracted automatically from sequence headers and used as identifiers in the results.
- Added the param `record_all_hits` to toggle logging every hit on or off. See `pipeline.run()` documentation for details.
