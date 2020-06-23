import random
from Bio import SeqIO
random.seed(43)


# Generates contigs with proteins of interest and CRISPR arrays
# This was used to generate the FASTA files in integration_data/end-to-end
# It does not need to be run or used anymore, it's here to document the process


def random_sequence(length: int) -> str:
    return "".join((random.choice("ACGT") for _ in range(length)))


def generate_crispr_array(repeat_type: str, repeat_count: int) -> str:
    repeat_sequence = "GTTTCAATCCACGCGCCCACGCGGGGCGCGAC" if repeat_type == 'cas12a' else "ATTTCAATCCACTCACCCATGAAGGGTGAGAC"
    spacer_length = random.choice(range(35, 42))
    spacers = [random_sequence(spacer_length) for _ in range(repeat_count)]
    return "".join([f"{repeat_sequence}{spacer}" for spacer in spacers])


def load_seqs():
    with open("integration_data/end-to-end/cas.fasta") as f:
        cas9, cas12a = [s.seq for s in SeqIO.parse(f, "fasta")]
    with open("integration_data/end-to-end/transposase.fasta") as f:
        transposase = list(SeqIO.parse(f, "fasta"))[0].seq
    data = {"nucleotides": (str(cas9), str(cas12a), str(transposase))}
    data["proteins"] = str(cas9[:-3].translate()), str(cas12a[:-3].translate()), str(transposase[:-3].translate())
    return data


def build_fasta_string(fasta, fastaid):
    return "\n".join([f">fasta{fastaid}_{n}\n{seq}" for n, seq in enumerate(fasta)])


if __name__ == '__main__':
    data = load_seqs()
    cas9, cas12a, transposase = data["nucleotides"]

    # fasta1
    # no features
    # cas9/transposase/CRISPR array
    # transposase and distant CRISPR array
    # cas9/CRISPR array
    # cas12/transposase/CRISPR array
    fasta1 = [random_sequence(8000),
              random_sequence(450) + cas9 + random_sequence(12) + transposase + random_sequence(31) + generate_crispr_array("cas9", 5) + random_sequence(7000),
              random_sequence(2301) + transposase + random_sequence(9000) + generate_crispr_array("cas12a", 10) + random_sequence(300),
              random_sequence(5) + cas9 + random_sequence(20) + generate_crispr_array("cas9", 2) + random_sequence(300),
              random_sequence(900) + cas12a + random_sequence(10) + transposase + random_sequence(10) + generate_crispr_array("cas12a", 20) + random_sequence(500)]
    # fasta2
    # CRISPR array
    # cas9
    fasta2 = [random_sequence(1700) + generate_crispr_array("cas9", 4) + random_sequence(2100),
              random_sequence(6703) + cas9 + random_sequence(3403)]
    # 3.fasta
    # 2 CRISPR arrays
    # no orfs
    # cas12a
    # no orfs
    # cas9/cas9/cas12a
    # no orfs
    fasta3 = [random_sequence(300) + generate_crispr_array("cas9", 7) + random_sequence(390) + generate_crispr_array("cas12a", 5) + random_sequence(800),
              random_sequence(4500),
              random_sequence(700) + generate_crispr_array("cas12a", 3) + random_sequence(9000),
              random_sequence(9000),
              random_sequence(500) + cas9 + random_sequence(30) + cas9 + random_sequence(35) + cas12a + random_sequence(500),
              random_sequence(90) + cas12a + random_sequence(10) + transposase + random_sequence(10) + generate_crispr_array("cas12a", 10) + random_sequence(400),
              random_sequence(8000)]

    for n, fasta in enumerate([fasta1, fasta2, fasta3]):
        version = n + 1
        output = build_fasta_string(fasta, version)
        with open(f"integration_data/end-to-end/fasta{version}.fasta", "w") as f:
            f.write(output)
