import csv
import itertools
import re
import io
from typing import List, Tuple, Optional
from Bio.Seq import Seq


class Feature(object):
    """ 
    Represents a gene or CRISPR repeat array. This is used internally by :class:`operon_analyzer.genes.Operon`, 
    but appears in the auto-generated documentation for reference. 
    """

    def __init__(self,
                 name: str,
                 coordinates: Tuple[int, int],
                 orfid: str,
                 strand: Optional[int],
                 accession: str,
                 e_val: Optional[float],
                 description: str,
                 sequence: str,
                 bit_score: Optional[int] = None,
                 raw_score: Optional[int] = None,
                 aln_len: Optional[int] = None,
                 pident: Optional[float] = None,
                 nident: Optional[float] = None,
                 mismatch: Optional[float] = None,
                 positive: Optional[float] = None,
                 gapopen: Optional[int] = None,
                 gaps: Optional[int] = None,
                 ppos: Optional[float] = None,
                 qcovhsp: Optional[int] = None):
        # Note: for CRISPR repeats, the pipeline does not identify the strand, and
        # since BLAST was not used, there is no e-value, so these values are set to
        # None.
        self.name = name
        self.coordinates = coordinates
        assert self.start < self.end
        self.orfid = orfid
        self.strand = strand
        self.accession = accession
        self.e_val = e_val
        self.bit_score = bit_score
        self.description = description
        self.sequence = sequence
        self.raw_score = raw_score
        self.aln_len = aln_len
        self.pident = pident
        self.nident = nident
        self.mismatch = mismatch
        self.positive = positive
        self.gapopen = gapopen
        self.gaps = gaps
        self.ppos = ppos
        self.qcovhsp = qcovhsp
        self.ignored_reasons = []

    @property
    def start(self) -> int:
        return self.coordinates[0]

    @property
    def end(self) -> int:
        return self.coordinates[1]

    def __hash__(self):
        return hash(f"{self.name}{self.start}{self.end}{self.sequence}")

    def __eq__(self, other: 'Feature'):
        return self.name == other.name and \
               self.start == other.start and \
               self.end == other.end and \
               self.orfid == other.orfid and \
               self.strand == other.strand and \
               self.accession == other.accession and \
               self.e_val == other.e_val and \
               self.bit_score == other.bit_score and \
               self.description == other.description and \
               self.sequence == other.sequence and \
               self.raw_score == other.raw_score and \
               self.aln_len == other.aln_len and \
               self.pident == other.pident and \
               self.nident == other.nident and \
               self.mismatch == other.mismatch and \
               self.positive == other.positive and \
               self.gapopen == other.gapopen and \
               self.gaps == other.gaps and \
               self.ppos == other.ppos and \
               self.qcovhsp == other.qcovhsp

    def ignore(self, reason: str):
        self.ignored_reasons.append(reason)

    def __len__(self):
        """ Length of feature in nucleotides """
        return self.end - self.start

    def __repr__(self) -> str:
        return f"<Feature {self.name} {self.start}..{self.end}>"


class Operon(object):
    """
    Provides access to features that were found in the same genomic region,
    which presumably comprise an actual operon. Whether this is true in reality
    must be determined by the user, if that is meaningful to them.
    """

    def __init__(self,
                 contig: str,
                 contig_filename: str,
                 start: int,
                 end: int,
                 features: List[Feature]):
        assert contig, "Missing contig name"
        assert 0 <= start and 0 <= end, "Invalid contig position"
        assert features, "Contig did not contain any features"
        self.contig = contig
        self.contig_filename = contig_filename
        self.start = start
        self.end = end
        self._features = features
        self._sequence = None

    def __hash__(self):
        feature_data = "".join([f"{feature.name}{feature.start}{feature.end}{feature.sequence}" for feature in self._features])
        return hash(f"{self.contig}{self.start}{self.end}{feature_data}")

    def __eq__(self, other):
        return self.contig == other.contig and \
               self.start == other.start and \
               self.end == other.end and \
               len(set(tuple(itertools.chain(self._features, other._features)))) == len(self._features)

    def set_sequence(self, sequence: Seq):
        """ Stores the nucleotide sequence of the operon. """
        self._sequence = sequence

    @property
    def feature_region(self) -> Tuple[int, int]:
        lower_bound, upper_bound = self.end, self.start
        for feature in self:
            lower_bound = min(feature.start, lower_bound)
            upper_bound = max(feature.end, upper_bound)
        return lower_bound, upper_bound

    @property
    def feature_region_sequence(self) -> str:
        """ Returns the nucleotide sequence of the operon, excluding the regions outside of the outermost Features. """
        lower_bound, upper_bound = self.feature_region
        return self._sequence[lower_bound:upper_bound]

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def all_genes(self):
        """ Iterates over all genes (i.e. not CRISPR arrays) in the operon regardless of whether it's been ignored. """
        yield from (feature for feature in self._features if feature.e_val is not None)

    @property
    def all_features(self):
        """ Iterates over all features in the operon regardless of whether it's been ignored. """
        yield from self._features

    def __iter__(self):
        """ Iterates over all features that haven't been judged to be unrelated to the operon. """
        yield from (feature for feature in self._features if not feature.ignored_reasons)

    def __len__(self):
        """ The number of relevant genes or CRISPR arrays in the operon. """
        return len(tuple(iter(self)))

    @property
    def feature_names(self):
        """ Iterates over the name of each feature in the operon """
        yield from (feature.name for feature in self)

    def get(self, feature_name: str, regex=False) -> List[Feature]:
        """ Returns a list of every Feature with a given name. """
        if regex:
            rx = re.compile(feature_name, re.IGNORECASE)
        else:
            rx = re.compile(f"^{feature_name}$", re.IGNORECASE)
        features = []
        for feature in self:
            if rx.search(feature.name):
                features.append(feature)
        return features

    def get_unique(self, feature_name: str, regex=False) -> Optional[Feature]:
        """ Returns a Feature or None if there is more than one Feature with the same name """
        features = self.get(feature_name, regex)
        if len(features) != 1:
            return None
        return features[0]

    def as_str(self) -> str:
        """
        Writes an Operon back out in the same CSV format that gene_finder produces. The text won't be
        completely identical in the case where floats have trailing decimals, or zero values in
        scientific format are recast as a simple float.
        """
        def optional_float(value, decimals: int) -> str:
            return f"{value:.{decimals}f}" if value else ""

        def optional_evalue(value) -> str:
            if value is None:
                return ""
            return f"{value:.2e}"

        def optional_float_or_int(value) -> str:
            if value is None:
                return ""
            if round(value) == value:
                return str(int(value))
            return f"{value:.1f}"

        def optional(value) -> str:
            return value if value not in ("", None) else ""

        output = io.StringIO(newline='')
        writer = csv.writer(output)
        for feature in self.all_features:
            feature_coords = f"{feature.start + 1}..{feature.end}" if feature.strand == 1 else f"{feature.end}..{feature.start + 1}"
            strand = feature.strand if feature.name != 'CRISPR array' else ''
            writer.writerow([f"{self.contig}",f"{self.start}..{self.end}",f"{feature.name}",f"{feature_coords}",f"{feature.orfid}",f"{optional(strand)}",f"{feature.accession}",f"{optional_evalue(feature.e_val)}",f"{feature.description}",f"{feature.sequence}",f"{optional_float_or_int(feature.bit_score)}",f"{optional(feature.raw_score)}",f"{optional(feature.aln_len)}",f"{optional_float(feature.pident, 3)}",f"{optional(feature.nident)}",f"{optional(feature.mismatch)}",f"{optional(feature.positive)}",f"{optional(feature.gapopen)}",f"{optional(feature.gaps)}",f"{optional_float(feature.ppos, 2)}",f"{optional(feature.qcovhsp)}",f"{self.contig_filename}"])
        return "%s" % output.getvalue().replace("\r", "")
