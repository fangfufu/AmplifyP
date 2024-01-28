# -*- coding: utf-8 -*-
"""Amplify P - Replication related."""

from typing import Dict, List
import logging

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

    def __init__(
        self,
        template: DNA,
        primer: Primer,
        settings: Settings,
    ) -> None:
        """Initialize the ReplicationConfig object."""
        self.padding_len = len(primer)

        self.primer = primer
        self.template = template

        self.template_seq: Dict[DNADirection, str] = {}
        # Add padding the 5' end of the template
        self.template_seq[DNADirection.FWD] = template.pad(self.padding_len).sequence
        # Add padding to the 3' end of the DNA, compute the complement.
        self.template_seq[DNADirection.REV] = (
            template.reverse().pad(self.padding_len).reverse().complement().sequence
        )

        logging.debug(
            f"Repliconf.__init__(): FWD: {self.template_seq[DNADirection.FWD]}"
        )
        logging.debug(
            f"Repliconf.__init__(): REV: {self.template_seq[DNADirection.REV]}"
        )

        self.settings = settings

        self.amplicon_start: List[int] = []
        self.amplicon_end: List[int] = []

    def clear(self) -> None:
        """Clears the lists of valid replication origins."""

        self.amplicon_start = []
        self.amplicon_end = []

    def range(self) -> range:
        """Return the valid search range for replication origin."""
        return range(
            0, len(self.template_seq[DNADirection.FWD]) - len(self.primer.sequence) + 1
        )

    def slice(self, i: int) -> slice:
        """Return a object for constructing ReplicationOrigin."""
        return slice(i, i + len(self.primer.sequence))

    def origin(self, direction: DNADirection, i: int) -> ReplicationOrigin:
        """Return the ith ReplicationOrigin."""
        return ReplicationOrigin(
            # We reverse the template sequence slice here, because Bill's
            # algorithm score from 3' to 5'.
            self.template_seq[direction][self.slice(i)][::-1],
            # We reverse the primer sequence here, if the template is 5'->3'
            # Otherwise the template is 3'->5', so we don't have to reverse
            # the primer
            self.primer.sequence[::-1] if direction else self.primer.sequence,
            self.settings,
        )

    def search(self) -> None:
        """Search for the valid replication origins in both directions."""
        self.clear()
        for direction in [DNADirection.FWD, DNADirection.REV]:
            logging.debug(f"Repliconf.search(): {direction}")
            for i in self.range():
                origin = self.origin(direction, i)
                if (
                    origin.primability() > self.settings.primability_cutoff
                    and origin.stability() > self.settings.stability_cutoff
                ):
                    logging.debug(f"Repliconf.search(): adding {origin}")
                    if direction:
                        self.amplicon_start.append(i - len(self.primer.sequence))
                    else:
                        self.amplicon_end.append(i + len(self.primer.sequence))

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
