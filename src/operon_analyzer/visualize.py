import os
import sys
from typing import Tuple, Dict, IO, List
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
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


def calculate_adjusted_operon_bounds(operon: Operon) -> Tuple[int, int]:
    """
    Calculates the offset that should be removed from each Feature's position
    as well as the total length of the operon.
    """
    low, high = sys.maxsize, 0
    # do one pass to find the bounds of the features we're interested in
    for feature in operon:
        low = min(low, feature.start)
        high = max(high, feature.end)
    assert low < high
    return low, high - low


def create_operon_figure(operon: Operon, feature_colors: dict):
    assert len(operon) > 0
    offset, operon_length = calculate_adjusted_operon_bounds(operon)

    graphic_features = []
    for feature in operon:
        if feature_colors is not None:
            graphic_feature = GraphicFeature(start=feature.start - offset,
                                            strand=feature.strand,
                                            end=feature.end - offset,
                                            label=feature.name,
                                            color=feature_colors[feature.name])
        else:
            graphic_feature = GraphicFeature(start=feature.start - offset,
                                            strand=feature.strand,
                                            end=feature.end - offset,
                                            label=feature.name)
        graphic_features.append(graphic_feature)

    record = GraphicRecord(sequence_length=operon_length,
                           features=graphic_features)

    ax, _ = record.plot(figure_width=5)
    record.plot(ax)
    return ax


def save_operon_figure(ax: Axes, out_filename: str):
    ax.figure.savefig(out_filename, bbox_inches='tight')
    plt.close()


def build_operon_dictionary(f: IO[str]) -> Dict[Tuple[str, int, int], Operon]:
    operons = {}
    lines = read_pipeline_output(f)
    for operon in assemble_operons(lines):
        operons[(operon.contig, operon.start, operon.end)] = operon
    return operons


def plot_operons(operons: List[Operon], output_directory: str, feature_colors: dict = None):
    for operon in operons:
        out_filename = build_image_filename(operon, output_directory)
        ax = create_operon_figure(operon, feature_colors)
        save_operon_figure(ax, out_filename)


if __name__ == '__main__':
    analysis_csv, pipeline_csv, image_directory = sys.argv[1:]
    good_operons = []

    with open(pipeline_csv) as f:
        operons = build_operon_dictionary(f)
    with open(analysis_csv) as f:
        for contig, start, end, result in load_analyzed_operons(f):
            if result != 'pass':
                continue
            op = operons.get((contig, start, end))
            if op is None:
                continue
            good_operons.append(op)
    plot_operons(good_operons, image_directory)
