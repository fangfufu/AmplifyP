# -*- coding: utf-8 -*-
"""Amplify P - Replication related."""

from typing import Dict, List

from .dna import DNA, Primer, DNADirection
from .origin import ReplicationOrigin
from .settings import Settings, DEFAULT_REPLICATION_CONFIG


class Repliconf:
    """
    Repliconf - a class representing a Replication Configuration.

    We define a "replication configuration" being a combination of a single
    primer and two padded target DNA sequence. An instance of this class should
    be created when we want switch the primer.

    The reason this class is needed is because primers vary by length, and the
    target sequence needs to be padded accordingly. There are two target
    sequences -- the forward sequence, and the reverse complement sequence.
    """

    def __init__(
        self,
        template: DNA,
        primer: Primer,
        settings: Settings = DEFAULT_REPLICATION_CONFIG,
    ) -> None:
        """Initialize the ReplicationConfig object."""
        self.padding_len = len(primer)

        self.primer = primer
        self.template = template

        # Bill's algorithm scores from 3' to 5'. We pre-reverse it here.
        self.primer_seq: str = primer.reverse().sequence

        # We don't pre-reverse these there in order to make things easier to
        # reason with.
        self.template_seq: Dict[DNADirection, str] = {}
        self.template_seq[DNADirection.FWD] = template.pad(self.padding_len).sequence
        self.template_seq[DNADirection.REV] = (
            template.reverse().complement().pad(self.padding_len).sequence
        )

        self.settings = settings

        self.amplicon_start: Dict[DNADirection, List[int]] = {}

    def clear(self) -> None:
        """
        Clears the lists of valid replication origins
        """

        self.amplicon_start[DNADirection.FWD] = []
        self.amplicon_start[DNADirection.REV] = []

    def range(self) -> range:
        """Return the valid search range for replication origin."""
        return range(
            0, len(self.template_seq[DNADirection.FWD]) - len(self.primer_seq) + 1
        )

    def slice(self, i: int) -> slice:
        """Return a slice of primer_seq for constructing ReplicationOrigin."""
        return slice(i, i + len(self.primer_seq))

    def origin(self, direction: DNADirection, i: int) -> ReplicationOrigin:
        """Return the ith ReplicationOrigin."""
        return ReplicationOrigin(
            # Note that we reverse the template sequence slice here, because
            # Bill's algorithm score from 3' to 5'.
            self.template_seq[direction][self.slice(i)][::-1],
            self.primer_seq,
            self.settings,
        )

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
                    self.amplicon_start[direction].append(
                        i if direction else i + len(self.primer_seq)
                    )

    def __eq__(self, other: object) -> bool:
        """
        Check if two ReplicationConfig objects are equal.

        Args:
            other (object): The object to compare with.

        Returns:
            bool: True if the objects are equal, False otherwise.
        """
        if not isinstance(other, Repliconf):
            return NotImplemented
        return self.primer == other.primer and self.template == other.template

    def __hash__(self) -> int:
        return hash((self.primer, self.template))

    def __str__(self) -> str:
        """Return the string representation of the ReplicationConfig object."""
        return f"ReplicationConfig: Primer: {self.primer}, Target: {self.template}"
