import pytest
import yaml, json, tempfile, os, shutil
from gene_finder.pipeline import Pipeline

def create_pipeline(conf_file):
    """Initialize a Pipeline object."""

    # load contents of yaml conf file
    stream = open(conf_file, 'r')
    conf = yaml.full_load(stream)
    p = Pipeline()
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
    return p

def merge_data():
    """
    Temporarily merge multiple sequences from
    the test data into a single fasta file.
    """

    # create a temporary directory for the merge file
    tmp = tempfile.TemporaryDirectory()
    merge = os.path.join(tmp.name, "merged_input.fasta")
    fasta1 = "tests/integration/integration_data/contigs/record_all_hits_test_1"
    fasta2 = "tests/integration/integration_data/contigs/record_all_hits_2"

    with open(merge,"wb") as f:
        with open(fasta1,"rb") as infile:
            shutil.copyfileobj(infile, f)
        with open(fasta2,"rb") as infile:
            shutil.copyfileobj(infile, f)
    
    return (tmp, merge)


@pytest.mark.slow
def test_with_blast():
    """Test pipeline using BLAST as the alignment tool.
    
    Use --runslow option to run from command line. 
    Otherwise, test will be skipped.
    """

    genomic_data = "tests/integration/integration_data/contigs/v_crass_J520_whole.fasta"
    conf_file = "tests/integration/configs/blast_integration_test.yaml"
    p = create_pipeline(conf_file)
    
    results = p.run(job_id="blast_test", data=genomic_data)
    
    assert len(results) == 1
    hits = results["NZ_CCKB01000071.1"]["Loc_78093-114093"]["Hits"]
    assert len(hits) == 11
    assert "Array_0" in hits
    os.remove("blast_test_results.csv")

@pytest.mark.slow
def test_multi_seq_fasta():
    """
    Test that the pipeline runs and produces the
    correct output for a multi-sequence fasta file.

    Use --runslow option to run from command line. 
    Otherwise, test will be skipped.
    """

    tmp, genomic_data = merge_data()
    conf_file = "tests/integration/configs/blast_integration_test.yaml"
    p = create_pipeline(conf_file)
    
    results = p.run(job_id="blast_test", data=genomic_data)
    
    assert len(results) == 2
    tmp.cleanup()
    os.remove("blast_test_results.csv")


@pytest.mark.slow
def test_gzip_fasta():
    """
    Test that the pipeline can open and correctly process
    a gzipped fasta file.
    """
    
    data = "tests/integration/integration_data/contigs/record_all_hits_test_1.gz"
    conf_file = "tests/integration/configs/blast_integration_test.yaml"
    p = create_pipeline(conf_file)
    results = p.run(job_id="blast_test", data=data, gzip=True)

    hits = results["KB405063.1"]["Loc_0-23815"]["Hits"]
    assert "Cas_all_hit-0" in hits
    assert "Array_0" in hits
    os.remove("blast_test_results.csv")
    

@pytest.mark.slow
def test_record_all_hits_1():
    """
    Testcase: input contig contains a putative
    CRISPR-Transposon system (tnsAB/tnsCD/cas genes/CRISPR array).
    
    Use --runslow option to run this test from the command line
    with pytest. 
    
    Otherwise, test will be skipped.
    """

    genomic_data = "tests/integration/integration_data/contigs/record_all_hits_test_1"
    conf_file = "tests/integration/configs/blast_integration_test.yaml"
    p = create_pipeline(conf_file)

    p.run("blast_test", data=genomic_data, record_all_hits=True)

    with open("blast_test_all_hits.json", "r") as f:
        hits = json.load(f)["KB405063.1"]
        assert len(hits) == 4
        assert "tnsAB" in hits
        assert "cas_all" in hits
        assert "tnsCD" in hits
        assert "CRISPR" in hits
    os.remove("blast_test_results.csv")
    os.remove("blast_test_all_hits.json")


@pytest.mark.slow
def test_record_all_hits_2():
    """
    Testcase: input contig has a hit for tnsA 
    but no hits for cas genes.
    """

    genomic_data = "tests/integration/integration_data/contigs/record_all_hits_2"
    conf_file = "tests/integration/configs/blast_integration_test.yaml"
    p = create_pipeline(conf_file)

    p.run("blast_test", data=genomic_data, record_all_hits=True)

    with open("blast_test_all_hits.json", "r") as f:
        hits = json.load(f)["JNYE01000019.1"]
        assert len(hits) == 1
        assert "tnsAB" in hits
        assert "cas_all" not in hits
        assert "tnsCD" not in hits
        assert "CRISPR" not in hits
    os.remove("blast_test_all_hits.json")
    os.remove("blast_test_results.csv")


@pytest.mark.slow
def test_pipeline_with_arbitrary_blast_flags1():
    """Tries to run the pipeline using arbitrary blast flags at each step."""

    genomic_data = "tests/integration/integration_data/contigs/record_all_hits_test_1"
    tnsAB_db = "tests/integration/integration_data/blast_databases/tnsAB/blast_db"
    cas_db = "tests/integration/integration_data/blast_databases/cas/blast_db"
    tns_CD = "tests/integration/integration_data/blast_databases/tnsCD/blast_db"
    p = Pipeline()
    p.add_seed_step(db=tnsAB_db, name="tnsAB", e_val=0.001, blast_type="PROT", word_size=4, gap_open=12)
    p.add_blast_step(db=tns_CD, name="tnsCD", e_val=0.001, blast_type="PSI", gap_trigger=24)
    p.add_filter_step(db=cas_db, name="cas", e_val=0.001, blast_type="PROT", lcase_masking=True)
    results = p.run("blast_test", data=genomic_data)
    # the test data should contain some genes of interest; check that we still found something
    assert len(results) != 0
    os.remove("blast_test_results.csv")


@pytest.mark.slow
def test_pipeline_with_arbitrary_blast_flags2():
    """Tries to run the pipeline using arbitrary blast flags at each step."""

    genomic_data = "tests/integration/integration_data/contigs/record_all_hits_test_1"
    tnsAB_db = "tests/integration/integration_data/blast_databases/tnsAB/blast_db"
    p = Pipeline()
    p.add_seed_with_coordinates_step(db=tnsAB_db, name="tnsAB", e_val=0.001, blast_type="PROT", word_size=4, gap_open=12)
    results = p.run("blast_test", data=genomic_data)
    # the test data should contain some genes of interest; check that we still found something
    assert len(results) != 0
    os.remove("blast_test_results.csv")
    