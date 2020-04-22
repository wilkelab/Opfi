import unittest
import tempfile, os, yaml, json
from crisposon.pipeline import Pipeline

class TestPipeline(unittest.TestCase):

    def test_with_blast(self):
        genomic_data = "tests/test_data/contigs/v_crass_J520_whole.fasta"
        conf_file = "tests/configs/blast_functional_test.yaml"

        # load contents of yaml conf file
        stream = open(conf_file, 'r')
        conf = yaml.load(stream)

        p = Pipeline(genome=genomic_data, id="v_crassostreae", min_prot_len=conf["min-prot-len"], span=conf["span"])
        for step in conf["steps"]:
            
            if step["type"] == "seed":
                p.add_seed_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"])
            
            elif step["type"] == "filter":
                p.add_filter_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"])
            
            elif step["type"] == "blast":
                p.add_blast_step(db=step["blast-db"], name=step["name"], e_val=step["e-val"], blast_type=step["blast-type"])
            
            else:
                p.add_crispr_step()
        
        results = p.run()
        
        self.assertEqual(len(results["Loc_78093-111502"]["Hits"]), 11)
        self.assertTrue("Array_0" in results["Loc_78093-111502"]["Hits"])
    