from collections import defaultdict
import itertools
import math
import os
import re
import sys
import random
import string
from typing import Tuple, Dict, IO, List, Optional, Iterable, Any
import matplotlib
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


def _get_feature_color(feature_name: str, feature_colors: Dict[str, Any]) -> Any:
    # Decides what color a feature should be. The keys in `feature_colors` are either strings
    # that should exactly match a feature name, or a regular expression. The values of the dictionary
    # should be anything that matplotlib interprets as a color.

    default_color = feature_colors.get("", "blue")
    # Since there was no exact match, either the key is a regular expression or the feature is one
    # the user doesn't care about, in which case we'll try to use the default
    for pattern, color in feature_colors.items():
        if re.search(pattern, feature_name, re.IGNORECASE):
            return color
    return default_color


def _find_colormap_bounds(operons: List[Operon], color_by_blast_statistic, other_operons: Optional[List[Operon]] = None) -> Tuple[float, float]:
    if other_operons is None:
        other_operons = []
    lower, upper = None, None
    if color_by_blast_statistic is not None:
        lower = sys.maxsize
        upper = -sys.maxsize
        for operon in itertools.chain(operons, other_operons):
            for feature in operon:
                if feature.strand is None:
                    # This is a CRISPR array
                    continue
                value = getattr(feature, color_by_blast_statistic)
                lower = min(value, lower)
                upper = max(value, upper)
    return lower, upper


def plot_operons(operons: List[Operon],
                 output_directory: str,
                 plot_ignored: bool = True,
                 color_by_blast_statistic: Optional[str] = None,
                 feature_colors: Optional[dict] = {},
                 nucl_per_line: Optional[int] = None,
                 show_accession: bool = False,
                 show_description: bool = False):
    """ 
    Takes :class:`operon_analyzer.genes.Operon` objects and saves plots of them to disk. 

    Args:
        operons (list): Operons to be plotted.
        output_directory (str): Path to the directory to save operon plots to.
        plot_ignored (bool, optional): Toggles plotting of features that were marked as ignorable 
            by :class:`operon_analyzer.rules.FilterSet` .
        color_by_blast_statistic (str, optional): Map an alignment quality statistic using the virdis
            color scale. For a list of alignment statistics captured by Opfi, see :ref:`opfi-output-format` .
        feature_colors (dict, optional): If a labeled database was used during candidate identification,
            features can be colored accordingly using "label": "feature-color" pairs. For more information
            about labeling sequence databases, see :ref:`labeling-sequences` .
        nucl_per_line (int, optional): Length (in base pairs) to wrap gene diagrams on. 
        show_accession (bool, optional): Show the accession number of the best hit for each plotted feature.
        show_description (bool, optional): Show the description of the best hit for each plotted feature.
    """
    lower, upper = _find_colormap_bounds(operons, color_by_blast_statistic)
    for operon in operons:
        out_filename = build_image_filename(operon, output_directory)
        fig = create_operon_figure(operon, plot_ignored, feature_colors, color_by_blast_statistic=color_by_blast_statistic, colormin=lower, colormax=upper, nucl_per_line=nucl_per_line, show_accession=show_accession, show_description=show_description)
        if fig is None:
            continue
        save_operon_figure(fig, out_filename)


def _calculate_paired_figure_dimensions(operon: Operon, other: Operon, operon_length: int, plot_ignored: bool):
    """ Determines the figure height and width needed to make a stacked operon figure.
    This is essentially a bunch of magic heuristics that seem to come up with values that
    minimize whitespace and avoid overlapping feature labels. """
    height_top = min(6, int(max(2, math.sqrt(len(operon)))))
    height_bottom = min(6, int(max(2, math.sqrt(len(other)))))
    figure_width = max(int(operon_length/1000), 8)
    return height_top, height_bottom, figure_width


def plot_operon_pairs(operons: List[Operon], other_operons: List[Operon], output_directory: str,
                      color_by_blast_statistic: Optional[str] = None,
                      plot_ignored: bool = False, feature_colors: Optional[dict] = {}):
    """ 
    Takes two lists of presumably related Operons, pairs them up such that the pairs overlap the same genomic region,
    and plots one on top of the other. This allows side-by-side comparison of two different pipeline runs, so that you can, for example,
    run your regular pipeline, then re-BLAST with a more general protein database like nr, and easily see how the annotations differ. 

    Args:
        operons (list): Operons to be plotted.
        other_operons (list): Related operons to be plotted for comparison.
        output_directory (str): Path to the directory to save operon plots to.
        plot_ignored (bool, optional): Toggles plotting of features that were marked as ignorable 
            by :class:`operon_analyzer.rules.FilterSet` .
        color_by_blast_statistic (str, optional): Map an alignment quality statistic using the virdis
            color scale. For a list of alignment statistics captured by Opfi, see :ref:`opfi-output-format` .
        feature_colors (dict, optional): If a labeled database was used during candidate identification,
            features can be colored accordingly using "label": "feature-color" pairs. For more information
            about labeling sequence databases, see :ref:`labeling-sequences` .
    """
    lower, upper = _find_colormap_bounds(operons, color_by_blast_statistic, other_operons=other_operons)
    for operon, other in make_operon_pairs(operons, other_operons):
        out_filename = build_image_filename(operon, output_directory)
        plot_operon_pair(operon, other, lower, upper, out_filename, color_by_blast_statistic, plot_ignored, feature_colors)


def plot_operon_pair(operon: Operon,
                     other: Operon,
                     lower: float,
                     upper: float,
                     out_filename: str,
                     color_by_blast_statistic: Optional[str] = None,
                     plot_ignored: bool = False,
                     feature_colors: Optional[dict] = {}):
    # Calculate the figure size and the range of coordinates in the contig that we will plot
    lower_coordinates_bound, operon_length = calculate_adjusted_operon_bounds(operon, plot_ignored)
    upper_coordinates_bound = lower_coordinates_bound + operon_length
    height_top, height_bottom, figure_width = _calculate_paired_figure_dimensions(operon, other, operon_length, plot_ignored)

    # We create a Figure and two Axes and make DNA Features Viewer
    # use them so that we can stay in control of the dimensions
    fig, (ax1, ax2) = plt.subplots(nrows=2,
                                   ncols=1,
                                   figsize=(figure_width, height_top + height_bottom),
                                   gridspec_kw={"height_ratios": (height_top, height_bottom)})

    # Plot the top operon
    create_operon_figure(operon,
                         plot_ignored,
                         feature_colors,
                         color_by_blast_statistic=color_by_blast_statistic,
                         colormin=lower,
                         colormax=upper,
                         existing_ax=ax1,
                         figure_height=height_top, show_description=True)

    # Plot the bottom operon. We set bounds here so that it exactly
    # matches the location in the contig of the top operon
    create_operon_figure(other,
                         plot_ignored,
                         feature_colors,
                         color_by_blast_statistic=color_by_blast_statistic,
                         colormin=lower,
                         colormax=upper,
                         bounds=(lower_coordinates_bound, upper_coordinates_bound),
                         existing_ax=ax2,
                         figure_height=height_bottom, show_description=True)

    # Save the figure to disk
    save_pair_figure(fig, out_filename)


def make_operon_pairs(operons: List[Operon], other: List[Operon]) -> List[Tuple[Operon, Operon]]:
    """ Takes two lists of operons and tries to find matching pairs that overlaps the same genomic region.
    The overlaps can be partial, or the span of one can be a subset of the other. The idea here is that
    if we run our pipeline on the same sequencing data with two different databases, the neighborhoods
    might not line up exactly, but we want to do direct comparisons of the same region (for example, you might
    run your seed/BLAST/filter database, and then want to re-BLAST with Swissprot or nr to refine the 
    matches).

    The regions covered by operons in `operons` are used as the reference point. """
    other_lookup = defaultdict(list)
    for operon in other:
        other_lookup[operon.contig].append(operon)

    pairs = []
    for operon in operons:
        candidates = other_lookup.get(operon.contig)
        if not candidates:
            continue
        best_overlap = 0
        best_candidate = None
        for candidate in candidates:
            overlap = _calculate_operon_overlap(operon, candidate)
            if overlap is None:
                continue
            if overlap >= 1.0:
                best_candidate = candidate
                best_overlap = overlap
                break
            if overlap > best_overlap:
                best_candidate = candidate
                best_overlap = overlap
        if best_candidate is not None:
            pairs.append((operon, best_candidate))
    return pairs


def _calculate_operon_overlap(operon: Operon, other_operon: Operon) -> Optional[float]:
    """ Calculates the fraction of an operon that overlaps with another operon. """
    operon_start, operon_end = min(operon.start, operon.end), max(operon.start, operon.end)
    other_operon_start, other_operon_end = min(other_operon.start, other_operon.end), max(other_operon.start, other_operon.end)

    if operon_start > other_operon_end or other_operon_start > operon_end:
        # these operons don't overlap at all
        return None
    # Find the lower of the two ends
    end = min(operon_end, other_operon_end)
    # Find the higher of the two starts
    start = max(operon_start, other_operon_start)
    # Determine how much overlap there is
    operon_length = operon_end - operon_start + 1
    return (end - start + 1) / operon_length


def create_operon_figure(operon: Operon,
                         plot_ignored: bool,
                         feature_colors: Optional[dict] = {},
                         bounds: Optional[Tuple[int, int]] = None,
                         existing_ax: Optional[Axes] = None,
                         color_by_blast_statistic: Optional[str] = None,
                         colormin: Optional[float] = None,
                         colormax: Optional[float] = None,
                         figure_height: Optional[int] = None,
                         nucl_per_line: Optional[int] = None,
                         show_accession: bool = False,
                         show_description: bool = False):
    """ Plots all the Features in an Operon. """
    if not plot_ignored and len(operon) == 0:
        return None

    if colormin is not None and colormax is not None:
        norm = matplotlib.colors.LogNorm(vmin=colormin + sys.float_info.min, vmax=colormax)
        cmap = matplotlib.cm.get_cmap('viridis_r')

    if not bounds:
        offset, operon_length = calculate_adjusted_operon_bounds(operon, plot_ignored)
    else:
        offset = bounds[0]
        operon_length = bounds[1] - bounds[0]

    graphic_features = []
    for feature in operon.all_features:
        if feature.ignored_reasons and not plot_ignored:
            continue

        if colormin is not None and colormax is not None and not feature_colors:
            value = getattr(feature, color_by_blast_statistic)
            if value is None:
                color = 'gray'
            elif value == 0:
                color = cmap(0)
            else:
                normed_value = norm(value)
                color = cmap(normed_value)
        else:
            color = _get_feature_color(feature.name, feature_colors)
        if bounds and (not bounds[0] <= feature.start or not bounds[1] >= feature.end):
            continue
        # we alter the name of CRISPR arrays to add the number of repeats
        # this is done here and not earlier in the pipeline so that it doesn't
        # affect any rules that need to match on the name
        name = feature.name
        if show_accession:
            name = f"{name} {feature.accession}"
        if show_description:
            name = f"{name} {feature.description}"
        strand = feature.strand
        if feature.name == "CRISPR array":
            copies, repeat, spacer = feature.description.split(",")
            _, count = copies.split()
            name = f"CRISPR array ({count})"
            strand = None  # don't let array have a directional arrow
        label = name if not feature.ignored_reasons else f"{name} (ignored)"
        graphic_feature = GraphicFeature(start=feature.start - offset,
                                         strand=strand,
                                         end=feature.end - offset,
                                         label=label,
                                         color=color)
        graphic_features.append(graphic_feature)
    record = GraphicRecord(sequence_length=operon_length,
                           features=graphic_features)

    if nucl_per_line is not None:
        fig, _ = record.plot_on_multiple_lines(nucl_per_line=nucl_per_line, max_label_length=255)
    else:
        figure_width = max(int(operon_length/900), 5)
        ax, _ = record.plot(figure_width=figure_width, figure_height=figure_height, ax=existing_ax, max_label_length=255)
        fig = ax.figure
    return fig


def save_operon_figure(fig, out_filename: str):
    """ Writes the operon figure to disk. """
    fig.savefig(out_filename, bbox_inches='tight')
    plt.close()


def save_pair_figure(fig, out_filename: str):
    """ Writes the operon pair figure to disk. """
    fig.savefig(out_filename, bbox_inches='tight')
    plt.close()


def build_operon_dictionary(f: IO[str]) -> Dict[Tuple[str, int, int], Operon]:
    """ Builds a dictionary of Operons since we need to be able to randomly access them. """
    operons = {}
    lines = read_pipeline_output(f)
    for operon in assemble_operons(lines):
        operons[(operon.contig, operon.contig_filename, operon.start, operon.end)] = operon
    return operons


def _load_passing_contigs(handle: IO):
    good_contigs = set()
    for contig, contig_filename, start, end, result in analyze.load_analyzed_operons(handle):
        if result[0] != 'pass':
            continue
        good_contigs.add((contig, contig_filename, start, end))
    return good_contigs


def make_clustered_stacked_operon_plots(operons: Iterable[Operon],
                                        other_operons: Iterable[Operon],
                                        image_directory: str,
                                        min_count: int = 10,
                                        plot_ignored: bool = False,
                                        color_by_blast_statistic: Optional[str] = None,
                                        feature_colors: Optional[dict] = None
                                        ):
    """
    Clusters operons and plots them on top of a reannotated version of the same operon. This allows the user to BLAST
    some set of data with a curated database, then re-BLAST it against a more general database, and compare the two
    directly in a cluster-specific manner.

    If a :class:`operon_analyzer.rules.FilterSet` was used during analysis, that same set should be evaluated on each operon before passing it
    into this function.

    Args:
        operons (iterable): The operons of interest.
        other_operons (iterable): Reannotated operons.
        image_directory (str): The directory where all subdirectories will be created. Will be created if it does not exist.
        min_count (int, optional): Groups must have at least this many systems in order to be plotted. Default is 10.
        plot_ignored (bool, optional): Toggles plotting of features that were marked as ignorable 
            by :class:`operon_analyzer.rules.FilterSet` .
        color_by_blast_statistic (str, optional): Map an alignment quality statistic using the virdis
            color scale. For a list of alignment statistics captured by Opfi, see :ref:`opfi-output-format` .
        feature_colors (dict, optional): If a labeled database was used during candidate identification,
            features can be colored accordingly using "label": "feature-color" pairs. For more information
            about labeling sequence databases, see :ref:`labeling-sequences` .
        nucl_per_line (int, optional): Length (in base pairs) to wrap gene diagrams on. 
    """
    if not os.path.exists(image_directory):
        os.makedirs(image_directory)

    clustered_operons = analyze.cluster_operons_by_feature_order(operons)
    _plot_clustered_stacked_operons(clustered_operons, other_operons, image_directory, plot_ignored=plot_ignored, color_by_blast_statistic=color_by_blast_statistic, feature_colors=feature_colors)


def make_clustered_operon_plots(analysis_csv: str,
                                operons: Iterable[Operon],
                                image_directory: str,
                                min_count: int = 10,
                                diff_against_csv: Optional[str] = None,
                                plot_ignored: bool = False,
                                color_by_blast_statistic: Optional[str] = None,
                                feature_colors: Optional[dict] = None,
                                nucl_per_line: Optional[int] = None
                                ):
    """ 
    Clusters operons by the order of their features and plots them in separate directories,
    adding the number of systems to the directory name. Only systems that passed the rules specified as a
    :class:`operon_analyzer.rules.RuleSet` object will be eligible to be plotted.

    If a :class:`operon_analyzer.rules.FilterSet` was used during analysis, that same FilterSet should be evaluated 
    on each operon before passing it into this function.

    Args:
        analysis_csv (str): Path to the CSV file created by :func:`operon_analyzer.analyze.analyze` .
        operons (iterable): The operons of interest.
        image_directory (str): The directory where all subdirectories will be created. Will be created if it does not exist.
        min_count (int, optional): Groups must have at least this many systems in order to be plotted. Default is 10.
        diff_against_csv (str): Path to a CSV file created by operon analyzer. Any clusters in this file 
            will be skipped when clustering operons from analysis_csv. The point of 
            this is to see only new systems when making slight alterations to rules.
        plot_ignored (bool, optional): Toggles plotting of features that were marked as ignorable 
            by :class:`operon_analyzer.rules.FilterSet` .
        color_by_blast_statistic (str, optional): Map an alignment quality statistic using the virdis
            color scale. For a list of alignment statistics captured by Opfi, see :ref:`opfi-output-format` .
        feature_colors (dict, optional): If a labeled database was used during candidate identification,
            features can be colored accordingly using "label": "feature-color" pairs. For more information
            about labeling sequence databases, see :ref:`labeling-sequences` .
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
    _plot_clustered_operons(plottable_operons, image_directory, plot_ignored, color_by_blast_statistic=color_by_blast_statistic, feature_colors=feature_colors, nucl_per_line=nucl_per_line)


def _plot_clustered_stacked_operons(clustered_operons: Dict[str, List[Operon]],
                                    other_operons: List[Operon],
                                    image_directory: str,
                                    plot_ignored: bool,
                                    color_by_blast_statistic: Optional[str] = None,
                                    feature_colors: Optional[dict] = None):
    for motif_items, ops in clustered_operons.items():
        motif_name = '-'.join(motif_items)
        motif_directory = _make_motif_directory_name(motif_name, len(ops), image_directory)
        os.makedirs(motif_directory)
        plot_operon_pairs(ops, other_operons, motif_directory, plot_ignored=plot_ignored, color_by_blast_statistic=color_by_blast_statistic, feature_colors=feature_colors)


def _plot_clustered_operons(clustered_operons: Dict[str, List[Operon]], image_dir: str, plot_ignored: bool, color_by_blast_statistic: Optional[str] = None, feature_colors: Optional[dict] = None, nucl_per_line: Optional[int] = None):
    """ Plots contigs, placing them in directories named by the count and motif of the operons.
    For example, if there are 527 operons with Cas9, glmS, a CRISPR array, and Cas1 (in that order or exactly reversed),
    the directory name will be 527-cas9-glms-array-cas1. """
    for motif_items, operons in clustered_operons.items():
        motif_name = '-'.join(motif_items)
        motif_directory = _make_motif_directory_name(motif_name, len(operons), image_dir)
        os.mkdir(motif_directory)
        plot_operons(operons, motif_directory, plot_ignored=plot_ignored, color_by_blast_statistic=color_by_blast_statistic, feature_colors=feature_colors, nucl_per_line=nucl_per_line)


def _make_motif_directory_name(motif_name: str, num_operons: int, image_dir: str) -> str:
    """ Makes a string to use as a directory name for some set of operons with the same motif """
    motif_name = motif_name.replace("CRISPR array", "array").replace(" ", "_")
    motif_name = f"{num_operons}-{motif_name}"
    # get the filename character limit for the file system we are trying to write to
    max_allowed_chars = os.statvfs(image_dir).f_namemax
    if len(motif_name) > max_allowed_chars:
        motif_name = motif_name[:max_allowed_chars - 10] + "-" + ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))
    return os.path.join(image_dir, motif_name)
