# -*- coding: utf-8 -*-
"""Amplify P - Replication related."""

from .dna import DNA, DNADirection
from .origin import ReplicationOrigin
from .settings import Settings


class ReplicationConfig:
    """
    A class representing a Replication Configuration.

    We define a "replication configuration" being a combination of a single
    primer and two padded target DNA sequence. An instance of this class should
    be created when we want switch the primer.

    The reason this class is needed is because primers vary by length, and the
    target sequence needs to be padded accordingly. There are two target
    sequences -- the forward sequence, and the reverse complement sequence.
    """

    def __init__(self, target: DNA, primer: DNA, settings: Settings) -> None:
        """Initialize the ReplicationConfig object."""
        min_overlap: int = settings.min_overlap
        padding_len: int = len(primer) - min_overlap

        self.primer = primer
        # We reverse the sequences here, because Bill's algorithm requires you
        # to count from the 3' end. This is documented in ReplicationOrigin's
        # class docstring.
        self.primer_sequence: str = primer.reverse().sequence
        self.forward_sequence: str = target.pad(padding_len).reverse().sequence
        self.reverse_sequence: str = (
            target.complement().pad(padding_len).reverse().sequence
        )

        self.settings = settings

    def range(self) -> range:
        """Return the range of the target in ReplicationConfig."""
        return range(0, len(self.forward_sequence) - len(self.primer_sequence) + 1)

    def slice(self, i: int) -> slice:
        """Return the slice of the target in ReplicationConfig."""
        return slice(i, i + len(self.primer_sequence))

    def origin(self, direction: DNADirection, i: int) -> ReplicationOrigin:
        """Return the origin of replication."""
        sequence = (
            self.forward_sequence
            if direction == DNADirection.FORWARD
            else self.reverse_sequence
        )
        return ReplicationOrigin(
            sequence[self.slice(i)],
            self.primer_sequence,
            self.settings,
        )
