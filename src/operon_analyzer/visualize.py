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
    """ Generates a name (and optionally a path) for the PNG we are going to produce. """
    filename = "{contig}-{start}-{end}.png".format(
                contig=operon.contig,
                start=operon.start,
                end=operon.end)
    if directory:
        return os.path.join(directory, filename)
    return filename


def calculate_adjusted_operon_bounds(operon: Operon, include_ignored: bool = True) -> Tuple[int, int]:
    """
    Calculates the offset that should be removed from each Feature's position
    as well as the total length of the operon.
    """
    low, high = sys.maxsize, 0
    # do one pass to find the bounds of the features we're interested in
    for feature in operon:
        if feature.ignored_reasons and not include_ignored:
            continue
        low = min(low, feature.start)
        high = max(high, feature.end)
    assert low < high
    return low, high - low


def create_operon_figure(operon: Operon, plot_ignored: bool):
    """ Plots all the Features in an Operon. """
    assert len(operon) > 0
    offset, operon_length = calculate_adjusted_operon_bounds(operon)
    graphic_features = []
    for feature in operon.all_features:
        if bool(feature.ignored_reasons) and not plot_ignored:
            continue
        label = feature.name if not feature.ignored_reasons else "{feature_name} (ignored)".format(feature_name=feature.name)
        graphic_feature = GraphicFeature(start=feature.start - offset,
                                         strand=feature.strand,
                                         end=feature.end - offset,
                                         label=label)
        graphic_features.append(graphic_feature)
    record = GraphicRecord(sequence_length=operon_length,
                           features=graphic_features)

    ax, _ = record.plot(figure_width=5)
    record.plot(ax)
    return ax


def save_operon_figure(ax: Axes, out_filename: str):
    """ Writes the operon figure to disk. """
    ax.figure.savefig(out_filename, bbox_inches='tight')
    plt.close()


def build_operon_dictionary(f: IO[str]) -> Dict[Tuple[str, int, int], Operon]:
    """ Builds a dictionary of Operons since we need to be able to randomly access them. """
    operons = {}
    lines = read_pipeline_output(f)
    for operon in assemble_operons(lines):
        operons[(operon.contig, operon.start, operon.end)] = operon
    return operons


def plot_operons(operons: List[Operon], output_directory: str, plot_ignored: bool = True):
    """ Takes Operons and saves plots of them to disk. """
    for operon in operons:
        out_filename = build_image_filename(operon, output_directory)
        ax = create_operon_figure(operon, plot_ignored)
        save_operon_figure(ax, out_filename)
