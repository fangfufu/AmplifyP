# -*- coding: utf-8 -*-
"""Amplify P - Replication related."""

from .dna import DNA
from .settings import Settings


class ReplicationConfig:
    """
    A class representing a Replication Configuration.

    We define a "replication configuration" being a combination of a single
    primer and two padded DNA target sequence. The reason this class is
    needed is because primers vary by length, and the target sequence needs to
    be padded accordingly. There are two target sequences, because we care about
    the forward sequence, and the reverse complement sequence.
    """

    def __init__(self, target: DNA, primer: DNA, config: Settings) -> None:
        """Initialize the ReplicationConfig object."""
        min_overlap: int = config.min_overlap
        padding_len: int = len(primer) - min_overlap

        self.primer = primer
        # We reverse the sequences here, because Bill's algorithm requires you
        # to count from the 3' end. This is documented in ReplicationOrigin's
        # class docstring.
        self.primer_sequence: str = target.reverse().sequence
        self.forward_sequence: str = target.pad(padding_len).reverse().sequence
        # We want the reverse of the "reverse complement", so we don't reverse.
        self.reverse_sequence: str = target.pad(padding_len).complement().sequence

    def range(self) -> range:
        """Return the range of the target in ReplicationConfig."""
        return range(0, len(self.forward_sequence) - len(self.primer_sequence) + 1)

    def slice(self, i: int) -> slice:
        """Return the slice of the target in ReplicationConfig."""
        return slice(i, i + len(self.primer_sequence))
