# -*- coding: utf-8 -*-
"""Amplify P - Replication related."""

from typing import Dict, List

from .dna import DNA, Primer, DNADirection
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

    def __init__(self, target: DNA, primer: Primer, settings: Settings) -> None:
        """Initialize the ReplicationConfig object."""
        min_overlap: int = settings.min_overlap
        padding_len: int = len(primer) - min_overlap

        self.primer = primer
        self.target = target

        # We reverse the sequences here, because Bill's algorithm requires you
        # to count from the 3' end. This is documented in ReplicationOrigin's
        # class docstring.
        self.target_seq: Dict[DNADirection, str] = {}
        self.target_seq[DNADirection.FWD] = target.pad(padding_len).reverse().sequence
        self.target_seq[DNADirection.REV] = (
            target.complement().pad(padding_len).reverse().sequence
        )
        self.primer_seq: str = primer.reverse().sequence

        self.settings = settings

    def __eq__(self, other: object) -> bool:
        """
        Check if two ReplicationConfig objects are equal.

        Args:
            other (object): The object to compare with.

        Returns:
            bool: True if the objects are equal, False otherwise.
        """
        if not isinstance(other, ReplicationConfig):
            return NotImplemented
        return self.primer == other.primer and self.target == other.target

    def __hash__(self) -> int:
        return hash((self.primer, self.target))

    def range(self) -> range:
        """Return the range of the target in ReplicationConfig."""
        return range(
            0, len(self.target_seq[DNADirection.FWD]) - len(self.primer_seq) + 1
        )

    def slice(self, i: int) -> slice:
        """Return the slice of the target in ReplicationConfig."""
        return slice(i, i + len(self.primer_seq))

    def origin(self, direction: DNADirection, i: int) -> ReplicationOrigin:
        """Return the origin of replication."""
        return ReplicationOrigin(
            self.target_seq[direction][self.slice(i)],
            self.primer_seq,
            self.settings,
        )

    def clear(self) -> None:
        """Clear the primer index."""
        for direction in DNADirection:
            self.primer.index.clear(self.target, direction)

    def search(self) -> None:
        """Search for the valid replication origins in both directions."""
        self.clear()
        for direction in DNADirection:
            for i in self.range():
                origin = self.origin(direction, i)
                if (
                    origin.primability() > self.settings.primability_cutoff
                    and origin.stability() > self.settings.stability_cutoff
                ):
                    self.primer.index.append(
                        self.target,
                        direction,
                        i,
                    )

    def __str__(self) -> str:
        """Return the string representation of the ReplicationConfig object."""
        return f"ReplicationConfig: Primer: {self.primer}, Target: {self.target}"


class PCRMixture:
    """
    A class representing a PCR Mixture.

    A PCRMixture represents the chemical environment of a PCR reaction. It
    contains one target DNA sequence, and a list of primers. PCRMixture
    generates ReplicationConfig.
    """

    def __init__(self, target: DNA, primers: List[Primer], settings: Settings) -> None:
        """Initialize the PCRMixture object."""
        self.target = target
        self.primers = primers
        self.settings = settings
        self.configs: List[ReplicationConfig] = []
        for primer in self.primers:
            self.configs.append(ReplicationConfig(self.target, primer, self.settings))

    def search(self) -> None:
        """Search for the valid replication origins in all configurations."""
        for config in self.configs:
            config.search()
