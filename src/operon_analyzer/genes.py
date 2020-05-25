from typing import List, Tuple, Optional


class Feature(object):
    """ Represents a gene or CRISPR repeat array. """

    def __init__(self,
                 name: str,
                 coordinates: Tuple[int, int],
                 orfid: str,
                 strand: Optional[int],
                 accession: str,
                 e_val: Optional[float],
                 description: str,
                 sequence: str):
        # Note: for CRISPR repeats, the pipeline does not identify the strand, and
        # since BLAST was not used, there is no e-value, so these values are set to
        # None.
        self.name = name
        self.coordinates = coordinates
        self.orfid = orfid
        self.strand = strand
        self.accession = accession
        self.e_val = e_val
        self.description = description
        self.sequence = sequence
        self.ignored_reasons = []

    def __hash__(self):
        return hash((self.name, self.coordinates, self.sequence))

    def __eq__(self, other: 'Feature') -> bool:
        return self.name == other.name and \
               self.coordinates == other.coordinates and \
               self.sequence == other.sequence

    def ignore(self, reason: str):
        self.ignored_reasons.append(reason)

    @property
    def start(self) -> int:
        """
        The leftmost position of the Feature, relative to the
        orientation of the Operon.
        """
        return min(self.coordinates)

    @property
    def end(self) -> int:
        """
        The rightmost position of the Feature, relative to the
        orientation of the Operon.
        """
        return max(self.coordinates)


class Operon(object):
    """
    Provides access to Features that were found in the same genomic region,
    which presumably comprise an actual operon. Whether this is true in reality
    must be determined by the user, if that is meaningful to them.
    """

    def __init__(self,
                 contig: str,
                 start: int,
                 end: int,
                 features: List[Feature]):
        assert contig, "Missing contig name"
        assert 0 <= start and 0 <= end, "Invalid contig position"
        assert features, "Contig did not contain any features"
        self.contig = contig
        self.start = start
        self.end = end
        self._features = features

    @property
    def all_features(self):
        yield from self._features

    def __iter__(self):
        yield from (feature for feature in self._features if not feature.ignored_reasons)

    def __len__(self):
        return len(self._features)

    @property
    def feature_names(self):
        """ Iterates over the name of each feature in the operon """
        yield from (feature.name for feature in self)

    def get(self, feature_name: str) -> List[Feature]:
        """ Returns a list of every Feature with a given name. """
        features = []
        for feature in self:
            if feature.name == feature_name:
                features.append(feature)
        return features

    def get_unique(self, feature_name: str) -> Optional[Feature]:
        """ Returns a Feature or None if there is more than one Feature with the same name """
        features = self.get(feature_name)
        if len(features) != 1:
            return None
        return features[0]
