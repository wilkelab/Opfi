import sys
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord


if __name__ == '__main__':
    sequencing_dir = sys.argv[1]
    figdir = sys.argv[2]

    accession, candidate_start, candidate_stop, contig, has_cas3, genecount = sys.stdin.read().strip().split()
    has_cas3 = False if has_cas3 == 'False' else True
    genecount = int(genecount)

    if genecount < 3:
        exit(0)

    candidate_start, candidate_stop = int(candidate_start), int(candidate_stop)
    features = []
    lower = sys.maxsize
    upper = 0
    try:
        with open("%s/classifications/%s/protein_position.csv" % (sequencing_dir, contig)) as f:
            genes = parse_position_file(f)
    except Exception as e:
        exit(0)

    for gene in genes:
        # find operon bounds
        lower = min(gene.ordered_start, gene.ordered_end, lower)
        upper = max(gene.ordered_start, gene.ordered_end, upper)

    if lower > upper:
        exit(0)
    for gene, color in zip(genes, colors):
        # if name == "CRISPR repeat":
        #     continue
        start = gene._start - lower
        end = gene._end - lower
        
        feat = GraphicFeature(start=start, strand=gene.strand, end=end, color=color, label=gene._name)
        features.append(feat)
    
    a = "-".join(accession.split("-")[:2])
    strand = 1 if candidate_start < candidate_stop else -1
    feat = GraphicFeature(start=candidate_start-lower, end=candidate_stop-lower, strand=strand, color="white", label=a)
    features.append(feat)
    if len(features) < 4:
        print("not enough features")
        exit(0)

    operon_length = int((upper - lower)*1.05)
    record = GraphicRecord(sequence_length=operon_length, features=features)

    ax, _ = record.plot(figure_width=5)
    record.plot(ax)
    ax.figure.savefig('{figdir}/{contig}-{acc}.png'.format(figdir=figdir, contig=contig, acc=a), bbox_inches='tight')
    plt.close()
