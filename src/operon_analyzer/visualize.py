import csv
import os
import sys
from typing import IO
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
from operon_analyzer.genes import Operon
from operon_analyzer.parse import parse_pipeline_results
from operon_analyzer.parse import _parse_coordinates


blue = (.365, .647, .855)
yellow = (.871, .812, .247)
green = (.376, .741, .408)
red = (.945, .345, .329)
gray = (.302, .302, .302)
orange = (.980, .643, .227)
pink = (.945, .486, .690)
brown = (.698, .569, .184)
purple = (.698, .463, .698)

colors = [blue, yellow, green, red, gray, orange, pink, brown, purple]


def load_analyzed_operons(f: IO):
    for line in csv.reader(filter(lambda line: not line.startswith("#"), f)):
        contig, coordinates, result = line
        start, end = _parse_coordinates(coordinates)
        yield contig, start, end, result


def generate_image_filename(operon: Operon, directory: str=None) -> str:
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

    for feature, color in zip(operon, colors):
        graphic_feature = GraphicFeature(start=feature.start - low,
                                         strand=feature.strand,
                                         end=feature.end - low,
                                         color=color,
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
    analysis_csv, pipeline_csv, image_directory = sys.argv[1:]
    operons = {}
    with open(pipeline_csv) as f:
        lines = csv.reader(f)
        for operon in parse_pipeline_results(lines):
            operons[(operon.contig, operon.start, operon.end)] = operon
    with open(analysis_csv) as f:
        for contig, start, end, result in load_analyzed_operons(f):
            if result != 'pass':
                continue
            operon = operons.get((contig, start, end))
            if operon is None:
                continue
            out_filename = generate_image_filename(operon, image_directory)
            create_operon_image(out_filename, operon)
