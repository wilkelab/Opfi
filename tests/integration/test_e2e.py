from io import StringIO
from crisposon.pipeline import Pipeline
from operon_analyzer.analyze import analyze
from operon_analyzer.rules import RuleSet, FilterSet


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


rs = RuleSet().require('cas12a')
output = StringIO()
analyze(pipeline_output, rs, output=output)
print(output.getvalue())
