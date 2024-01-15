# -*- coding: utf-8 -*-
"""Amplify P - Replication related."""

from typing import Dict, List

from .dna import DNA, Primer, DNADirection
from .origin import ReplicationOrigin
from .settings import Settings


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

    def __init__(self, template: DNA, primer: Primer, settings: Settings) -> None:
        """Initialize the ReplicationConfig object."""
        self.padding_len = len(primer)

        self.primer = primer
        self.template = template

        # We reverse the sequences here, because Bill's algorithm requires you
        # to count from the 3' end. This is documented in ReplicationOrigin's
        # class docstring.
        self.template_seq: Dict[DNADirection, str] = {}
        self.template_seq[DNADirection.FWD] = (
            template.pad(self.padding_len).reverse().sequence
        )
        self.template_seq[DNADirection.REV] = (
            template.complement().pad(self.padding_len).reverse().sequence
        )
        self.primer_seq: str = primer.reverse().sequence

        self.settings = settings

        self.origin_start: Dict[DNADirection, List[int]] = {}

    def clear_origin_start(self) -> None:
        """
        Clears the valid origin dictionary.

        The valid origin dictionary stores the valid origins for each DNA direction.
        After calling this method, the dictionary will be empty.
        """

        self.origin_start[DNADirection.FWD] = []
        self.origin_start[DNADirection.REV] = []

    def idx_repliconf_to_template(self, repliconf_idx: int) -> int:
        """
        Repliconf / Template index conversion.

        Converts an index that is relative the repliconf DNA sequence (which is
        reverse and padded) to an index relative to the template sequence.

        Args:
            origin_idx (int): The index in the repliconf sequence.

        Returns:
            int: The corresponding index in the template sequence.
        """
        repliconf_idx -= self.padding_len
        return len(self.template) - repliconf_idx - len(self.primer)

    def range(self) -> range:
        """Return the valid origin range for this ReplicationConfig."""
        return range(
            0, len(self.template_seq[DNADirection.FWD]) - len(self.primer_seq) + 1
        )

    def slice(self, i: int) -> slice:
        """Return a slice of for constructing ReplicationOrigin."""
        return slice(i, i + len(self.primer_seq))

    def origin(self, direction: DNADirection, i: int) -> ReplicationOrigin:
        """Return an origin of replication."""
        return ReplicationOrigin(
            self.template_seq[direction][self.slice(i)],
            self.primer_seq,
            self.settings,
        )

    def search(self) -> None:
        """Search for the valid replication origins in both directions."""
        # self.clear()
        for direction in DNADirection:
            for i in self.range():
                origin = self.origin(direction, i)
                if (
                    origin.primability() > self.settings.primability_cutoff
                    and origin.stability() > self.settings.stability_cutoff
                ):
                    self.origin_start[direction].append(
                        self.idx_repliconf_to_template(i)
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
