import csv
import os
import sys
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
from operon_analyzer.analyze import load_analyzed_operons
from operon_analyzer.genes import Operon
from operon_analyzer.parse import assemble_operons, read_pipeline_output


def build_image_filename(operon: Operon, directory: str = None) -> str:
    filename = "{contig}-{start}-{end}.png".format(
                contig=operon.contig,
                start=operon.start,
                end=operon.end)
    if directory:
        return os.path.join(directory, filename)
    return filename


def create_operon_image(out_filename: str, operon: Operon):
    graphic_features = []
    low, high = sys.maxsize, 0
    # do one pass to find the bounds of the features we're interested in
    for feature in operon:
        low = min(low, feature.start)
        high = max(high, feature.end)

    for feature in operon:
        graphic_feature = GraphicFeature(start=feature.start - low,
                                         strand=feature.strand,
                                         end=feature.end - low,
                                         label=feature.name)
        graphic_features.append(graphic_feature)

    if not graphic_features:
        return False

    record = GraphicRecord(sequence_length=high-low,
                           features=graphic_features)

    ax, _ = record.plot(figure_width=5)
    record.plot(ax)
    ax.figure.savefig(out_filename, bbox_inches='tight')
    plt.close()
    return True


if __name__ == '__main__':
    # TODO: This makes too many assumptions about what the user wants
    # TODO: what if the user wants to plot failed contigs?
    # TODO: split this functionality out and make it easier to choose what to plot
    analysis_csv, pipeline_csv, image_directory = sys.argv[1:]
    operons = {}
    with open(pipeline_csv) as f:
        lines = read_pipeline_output(f)
        for operon in assemble_operons(lines):
            operons[(operon.contig, operon.start, operon.end)] = operon
    with open(analysis_csv) as f:
        for contig, start, end, result in load_analyzed_operons(f):
            if result != 'pass':
                continue
            op = operons.get((contig, start, end))
            if op is None:
                continue
            out_filename = build_image_filename(op, image_directory)
            create_operon_image(out_filename, op)
