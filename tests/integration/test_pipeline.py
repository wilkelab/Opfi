import pytest
import yaml
from crisposon.pipeline import Pipeline

@pytest.mark.slow
def test_with_blast():
    """Test pipeline using BLAST as the alignment tool.
    
    Use --runslow option to run from command line. 
    Otherwise, test will be skipped.
    """

    genomic_data = "tests/integration/integration_data/contigs/v_crass_J520_whole.fasta"
    conf_file = "tests/integration/configs/blast_integration_test.yaml"

    # load contents of yaml conf file
    stream = open(conf_file, 'r')
    conf = yaml.load(stream)

    p = Pipeline(genome=genomic_data, id="v_crassostreae", 
                min_prot_len=conf["min-prot-len"], span=conf["span"])
    for step in conf["steps"]:
        
        if step["type"] == "seed":
            p.add_seed_step(db=step["blast-db"], name=step["name"], 
                                e_val=step["e-val"], blast_type=step["blast-type"])
        
        elif step["type"] == "filter":
            p.add_filter_step(db=step["blast-db"], name=step["name"],
                                e_val=step["e-val"], blast_type=step["blast-type"])
        
        elif step["type"] == "blast":
            p.add_blast_step(db=step["blast-db"], name=step["name"], 
                                e_val=step["e-val"], blast_type=step["blast-type"])
        
        else:
            p.add_crispr_step()
    
    results = p.run()
    
    assert len(results["Loc_78093-111502"]["Hits"]) == 11
    