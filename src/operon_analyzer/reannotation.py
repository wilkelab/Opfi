import re
from collections import defaultdict
from typing import DefaultDict, Dict, List, Optional, Tuple

from operon_analyzer import analyze, genes

ReannotationCounts = Dict[Tuple[str, int], Dict[str, int]]
ReannotationFractions = Dict[Tuple[str, int], Dict[str, float]]


def summarize(operons: List[genes.Operon], reannotated_operons: List[genes.Operon]):
    """ For each protein in each operon, prints out how frequently it was reannotated to a particular protein.
    For example, a putative Cas9 may be identified as a transposase when BLASTed with a more expansive database. """
    clusters, reblasted_operons = _prepare_operons_for_counting(operons, reannotated_operons)
    for label, cloperons in clusters.items():
        counts = _count_cluster_reannotations(cloperons, reblasted_operons)
        if not counts:
            continue
        fractions = _convert_reannotation_counts_to_fractions(counts)
        for (feature, count), reannotations in fractions.items():
            print(_format_reannotation_summary("-".join(label), len(cloperons), fractions))


def _count_cluster_reannotations(operons: List[genes.Operon],
                                 reannotated_operons: Dict[str, genes.Operon]) -> ReannotationCounts:
    """
    After clustering operons by motif, the next question we'd often like to ask is how good
    the annotations are. This function takes the operons from one cluster, and operons
    that have been re-BLASTed against a more general database, and determines how often
    the annotations from the first round of BLASTing were reannotated in a parcticular way
    during the second round.

    For example, in a cluster with 130 operons, a particular region of DNA is annotated as cas9 using a curated
    database, and after BLASTing against nr, we find:

    - 107 S. pyogenes Cas9
    - 13 S. pyogenes DNA methyltransferase
    - 8 Unknown hypothetical protein
    - 2 Chloride intracellular channel protein 1

    While a small number of annotations look like other proteins, they're probably just noise, and we may
    tentatively conclude that this is a bonafide Cas9.

    After calculating the counts, this function then converts them to rough fractions:
    - 82 S. pyogenes Cas9
    - 10 S. pyogenes DNA methyltransferase
    - 6 Unknown hypothetical protein
    - 1 Chloride intracellular channel protein 1

    Partial overlaps with fraction greater than min_overlap_threshold coverage of the original operon
    will be included in the counts.

    """
    cluster_reannotations = defaultdict(dict)

    for operon in operons:
        reannotated_operon = reannotated_operons.get(operon.contig)
        if reannotated_operon is None:
            # We don't have a reannotation for this particular operon
            # This is common if we only re-BLAST a subset of our data
            # for efficiency's sake
            continue

        # Count the reannotations for a single pair of operons and merge those counts into the 
        # running total for the entire cluster
        update = _count_reannotations(operon, reannotated_operon)
        for feature_name, reannotation_data in update.items():
            for reannotated_feature_name, count in reannotation_data.items():
                current_count = cluster_reannotations[feature_name].get(reannotated_feature_name, 0)
                cluster_reannotations[feature_name][reannotated_feature_name] = count + current_count

            # If we get no reannotations for feature_name, we still need to add feature_name to the
            # cluster_reannotations dictionary, so that we can see when features are NEVER reannotated
            if not cluster_reannotations.get(feature_name):
                cluster_reannotations[feature_name] = defaultdict(int)
    return cluster_reannotations


def _count_reannotations(operon: genes.Operon,
                        reannotated_operon: genes.Operon):
    # To disambiguate features that occur more than once in each operon, we keep track of how many we have
    # of each. This way, we can refer to "transposase" and "transposase-2" as separate entities that
    # might BLAST to two completely different proteins.
    reannotations = {}
    feature_counts = defaultdict(int)
    for feature in operon:
        # determine which feature we're referring to in case more than one has the same name
        # if there is more than one, figure out what the number we need to append to the name should be
        count = feature_counts.get(feature.name)
        if not count:
            feature_name = (feature.name, 1)
            reannotations[feature_name] = defaultdict(int)
            feature_counts[feature.name] = 1
        else:
            feature_counts[feature.name] += 1
            feature_name = (feature.name, count + 1)
            if feature_name not in reannotations:
                reannotations[feature_name] = defaultdict(int)

        for reannotated_feature in reannotated_operon:
            if feature.start <= reannotated_feature.start and feature.end >= reannotated_feature.end:
                # reannotated_feature is completely overlapped by feature. It's possible
                # that there are multiple reannotated hits that map to this feature.
                # if that is the case, we could get >100% identities
                reannotations[feature_name][reannotated_feature.name] += 1
            elif reannotated_feature.start <= feature.start and reannotated_feature.end >= feature.end:
                # feature is completely overlapped by reannotated_feature
                reannotations[feature_name][reannotated_feature.name] += 1
    return reannotations


def _convert_reannotation_counts_to_fractions(reannotations: ReannotationCounts) -> ReannotationFractions:
    """ Normalizes raw counts to fractions of all reannotations. """
    fractional_reannotations = {}
    for feature_name, count_data in reannotations.items():
        fractional_reannotations[feature_name] = {}
        if not count_data:
            continue
        else:
            total = sum(count_data.values())
            for other_name, count in sorted(count_data.items(), key=lambda x: -x[1]):
                fraction = count / total
                fractional_reannotations[feature_name][other_name] = fraction
    return fractional_reannotations


def _format_reannotation_summary(contig_label: str,
                                 num_operons: int,
                                 fractional_reannotations: ReannotationCounts,
                                 show_top: Optional[int] = 3) -> str:
    """ Converts summary data into a formatted string that can be printed. """
    output = [f"{num_operons}-{contig_label}"]
    for (feature_name, feature_count), reannotated_feature_data in fractional_reannotations.items():
        formatted_feature_name = feature_name if feature_count == 1 else f"{feature_name}-{feature_count}"
        output.append(f"  {formatted_feature_name}:")
        limit = len(reannotated_feature_data) if show_top is None else show_top
        for _, (reannotated_feature_name, fraction) in zip(range(limit), reannotated_feature_data.items()):
            percentage = 100 * fraction
            output.append(f"    {percentage :4.0f}  {reannotated_feature_name}")
    return "\n".join(output)


def _prepare_operons_for_counting(operons: List[genes.Operon], compare_to_operons: List[genes.Operon]) \
        -> Tuple[DefaultDict[Tuple[str, ...], List[genes.Operon]], Dict[str, genes.Operon]]:
    """ Takes two lists of operons and reformats them so that we can perform our reannotation analysis. """
    clusters = analyze.cluster_operons_by_feature_order(operons)
    reannotated_operons = {operon.contig: operon for operon in compare_to_operons}
    return clusters, reannotated_operons


def _has_at_least_one_feature_with_fractional_matches(summary: ReannotationFractions,
                                                     feature_name: str,
                                                     regex_pattern: str,
                                                     min_matching_threshold: float):
    """
    Determines if the cluster has a feature that is reannotated to something that
    matches regex_pattern more frequently than min_matching_threshold.
    """
    assert 0 < min_matching_threshold <= 1.0
    rx = re.compile(regex_pattern, re.IGNORECASE)
    for (feature_name, _), fractions in summary.items():
        total = 0
        for other_feature_name, fraction in fractions.items():
            if rx.search(other_feature_name):
                total += fraction
        if total > min_matching_threshold:
            return True
    return False
