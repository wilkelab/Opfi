import os
import sys
from typing import Tuple, Dict, IO, List, Optional, Iterable, Set
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from dna_features_viewer import GraphicFeature, GraphicRecord
from operon_analyzer.genes import Operon
from operon_analyzer.parse import assemble_operons, read_pipeline_output
from operon_analyzer import analyze


ContigDescriptor = Tuple[str, str, int, int]  # contig accession ID, contig filename, start coordinate, end coordinate


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
    for feature in operon.all_features:
        if feature.ignored_reasons and not include_ignored:
            continue
        low = min(low, feature.start)
        high = max(high, feature.end)
    assert low < high
    return low, high - low


def create_operon_figure(operon: Operon, plot_ignored: bool, feature_colors: Optional[dict] = {}):
    """ Plots all the Features in an Operon. """
    if not plot_ignored and len(operon) == 0:
        return None
    # set the default color to the user-supplied one, if given, otherwise use blue
    default_color = feature_colors.get("", "blue")

    offset, operon_length = calculate_adjusted_operon_bounds(operon, plot_ignored)
    graphic_features = []
    for feature in operon.all_features:
        if feature.ignored_reasons and not plot_ignored:
            continue
        color = feature_colors.get(feature.name, default_color)
        # we alter the name of CRISPR arrays to add the number of repeats
        # this is done here and not earlier in the pipeline so that it doesn't
        # affect any rules that need to match on the name
        if feature.name == "CRISPR array":
            copies, repeat, spacer = feature.description.split(",")
            _, count = copies.split()
            feature.name = f"CRISPR array ({count})"
        label = feature.name if not feature.ignored_reasons else f"{feature.name} (ignored)"
        graphic_feature = GraphicFeature(start=feature.start - offset,
                                         strand=feature.strand,
                                         end=feature.end - offset,
                                         label=label,
                                         color=color)
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
        operons[(operon.contig, operon.contig_filename, operon.start, operon.end)] = operon
    return operons


def plot_operons(operons: List[Operon], output_directory: str, plot_ignored: bool = True, feature_colors: Optional[dict] = {}):
    """ Takes Operons and saves plots of them to disk. """
    for operon in operons:
        out_filename = build_image_filename(operon, output_directory)
        ax = create_operon_figure(operon, plot_ignored, feature_colors)
        if ax is None:
            continue
        save_operon_figure(ax, out_filename)


def _load_passing_contigs(handle: IO):
    good_contigs = set()
    for contig, contig_filename, start, end, result in analyze.load_analyzed_operons(handle):
        if result[0] != 'pass':
            continue
        good_contigs.add((contig, contig_filename, start, end))
    return good_contigs


def make_clustered_operon_plots(analysis_csv: str,
                                operons: Iterable[Operon],
                                image_directory: str,
                                min_count: int = 10,
                                diff_against_csv: Optional[str] = None,
                                plot_ignored: bool = False,
                                feature_colors: Optional[dict] = None
                                ):
    """ Clusters operons by the order of their features and plots them in separate directories,
    adding the number of systems to the directory name. Only systems that passed the Ruleset will
    be eligible to be plotted.

    If a FilterSet was used during analysis, that same FilterSet should be evaluated on each operon before passing it
    into this function.

    analysis_csv:       path to the CSV file created by operon analyzer
    operons:            the operons of interest
    image_directory:    the directory where all subdirectories will be created. Will be created if it does not exist
    min_count:          groups must have at least this many systems in order to be plotted
    diff_against_csv:   path to a CSV file created by operon analyzer. Any clusters in this file
                        will be skipped when clustering operons from analysis_csv. The point of
                        this is to see only new systems when making slight alterations to rules
    plot_ignored:       plot ignored features
    feature_colors:     a dictionary of feature names and their colors
    """
    with open(analysis_csv) as f:
        passing_contigs = _load_passing_contigs(f)

    diff_contigs = set()
    if diff_against_csv is not None:
        with open(diff_against_csv) as f:
            diff_contigs = _load_passing_contigs(f)

    if not os.path.exists(image_directory):
        os.mkdir(image_directory)

    plottable_operons = []
    diff_operons = []

    for operon in operons:
        key = operon.contig, operon.contig_filename, operon.start, operon.end
        if key in passing_contigs:
            plottable_operons.append(operon)
        if key in diff_contigs:
            diff_operons.append(operon)

    all_clustered_operons = analyze.cluster_operons_by_feature_order(plottable_operons)
    diff_clustered_operon_motifs = set(analyze.cluster_operons_by_feature_order(diff_operons).keys())
    plottable_operons = {k: ops for k, ops in all_clustered_operons.items()
                         if k not in diff_clustered_operon_motifs and len(ops) >= min_count}
    _plot_clustered_operons(plottable_operons, image_directory, plot_ignored, feature_colors)


def _plot_clustered_operons(clustered_operons: Dict[str, List[Operon]], image_dir: str, plot_ignored: bool, feature_colors: Optional[dict]):
    """ Plots contigs, placing them in directories named by the count and motif of the operons.
    For example, if there are 527 operons with Cas9, glmS, a CRISPR array, and Cas1 (in that order or exactly reversed),
    the directory name will be 527-cas9-glms-array-cas1. """
    for motif_items, operons in clustered_operons.items():
        motif_name = '-'.join(motif_items)
        motif_directory = _make_motif_directory_name(motif_name, len(operons), image_dir)
        os.mkdir(motif_directory)
        plot_operons(operons, motif_directory, plot_ignored=plot_ignored, feature_colors=feature_colors)


def _make_motif_directory_name(motif_name: str, num_operons: int, image_dir: str) -> str:
    """ Makes a string to use as a directory name for some set of operons with the same motif """
    motif_name = motif_name.replace("CRISPR array", "array").replace(" ", "_")
    motif_name = f"{num_operons}-{motif_name}"
    return os.path.join(image_dir, motif_name)
