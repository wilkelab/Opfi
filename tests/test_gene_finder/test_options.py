from gene_finder.option_handling import build_blastp_command, build_psiblast_command

def test_build_command1():
    # simply tests that the subprocess argument list is being constructed as expected
    kwargs = {"num_threads": 4, "ungapped": True, "qcov_hsp_perc": 90}
    cmd = build_blastp_command(query="my_query", db="my_db", evalue=0.001, kwargs=kwargs, 
                               output_fields="qseqid sseqid stitle", out="out.tsv", blastp_path='blastp')
    assert cmd == ["blastp", "-query", "my_query", "-db", "my_db", "-evalue", "0.001", "-out", "out.tsv",
                   "-num_threads", "4", "-ungapped", "-qcov_hsp_perc", "90", "-outfmt", "6 qseqid sseqid stitle"]


def test_build_command2():
    # tests that unknown/prohibited arguments are ignored
    kwargs = {"not_a_real": "blast_argument", "outfmt": 2}
    cmd = build_blastp_command(query="my_query", db="my_db", evalue=0.001, kwargs=kwargs, 
                               output_fields="qseqid sseqid stitle", out="out.tsv", blastp_path='blastp')
    assert cmd == ["blastp", "-query", "my_query", "-db", "my_db", "-evalue", "0.001", "-out", "out.tsv",
                   "-outfmt", "6 qseqid sseqid stitle"]


def test_build_command3():
    # test the psi blast version of this, for good measure
    kwargs = {"num_iterations": 3, "gap_trigger": 24}
    cmd = build_psiblast_command(query="my_query", db="my_db", evalue=0.001, kwargs=kwargs, 
                               output_fields="qseqid sseqid stitle", out="out.tsv", psiblast_path='psiblast')
    assert cmd == ["psiblast", "-query", "my_query", "-db", "my_db", "-evalue", "0.001", "-out", "out.tsv",
                   "-num_iterations", "3", "-gap_trigger", "24", "-outfmt", "6 qseqid sseqid stitle"]
