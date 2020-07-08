import os
import shutil
import tempfile
import pytest
from io import StringIO
from gene_finder.pipeline import Pipeline
from operon_analyzer.rules import RuleSet, FilterSet
from operon_analyzer.analyze import analyze, load_analyzed_operons
from operon_analyzer.visualize import build_operon_dictionary, plot_operons


def make_pngs(text, operons):
    tempdir = tempfile.mkdtemp()
    try:
        good_operons = []
        for contig, contig_filename, start, end, result in load_analyzed_operons(text.strip().split('\n')):
            if result != ['pass']:
                continue
            op = operons.get((contig, contig_filename, start, end))
            if op is None:
                continue
            good_operons.append(op)
        plot_operons(good_operons, tempdir)
        files = os.listdir(tempdir)
        count = len([f for f in files if f.endswith(".png")])
        return count
    except Exception as e:
        raise e
    finally:
        # clean up the directory and any PNGs that were made
        shutil.rmtree(tempdir)


@pytest.mark.slow
@pytest.mark.parametrize('requirement,pass_count,fail_count', [
    ('cas12a', 2, 1),
    ('cas9', 1, 2),
    ])
def test_end_to_end(requirement, pass_count, fail_count):
    pipeline_output = []

    for n in range(1, 4):
        filename = f"tests/integration/integration_data/end-to-end/fasta{n}.fasta"
        p = Pipeline()
        p.add_seed_step(db="tests/integration/integration_data/end-to-end/cas.prot", name="cas", e_val=0.000001, blast_type="PROT")
        p.add_filter_step(db="tests/integration/integration_data/end-to-end/transposase.prot", name="transposases", e_val=0.000001, blast_type="PROT")
        p.add_crispr_step()

        outfile = f"/tmp/e2e{n}.csv"
        results = p.run(filename, outfrmt="CSV", outfile=outfile, min_prot_len=30, span=10000)

        with open(outfile) as f:
            for line in f.readlines():
                pipeline_output.append(line)

    operons = build_operon_dictionary(pipeline_output)

    rs = RuleSet().require(requirement)
    fs = FilterSet().must_be_within_n_bp_of_anything(1000)
    output = StringIO()
    analyze(pipeline_output, rs, fs, output=output)

    text = output.getvalue()
    assert text.count('pass') == pass_count
    assert text.count('fail') == fail_count
    image_count = make_pngs(text, operons)
    assert image_count == pass_count
