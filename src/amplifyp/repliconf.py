"""Replication configuration-related classes for AmplifyP."""

import logging
from dataclasses import dataclass
from functools import cached_property

from .dna import DNA, DNADirection, Primer
from .origin import ReplicationOrigin
from .settings import Settings


@dataclass(slots=True)
class OriginIdx:
    """A class for storing the locations of replication origins.

    This class holds the indices of forward and reverse replication origins
    found within a DNA template. It also tracks whether a search has been
    performed.

    Attributes:
        fwd (List[int]): A list of indices for forward replication origins.
        rev (List[int]): A list of indices for reverse replication origins.
        searched (bool): A flag indicating whether a search for origins has
                         been completed.
    """

    fwd: list[int]
    rev: list[int]
    searched: bool

    def clear(self) -> None:
        """Clear the origin index, resetting all lists and flags."""
        self.fwd = []
        self.rev = []
        self.searched = False


class Repliconf:
    """A class representing a replication configuration.

    A "replication configuration" is defined as a combination of a single
    primer and a padded target DNA sequence. This class is instantiated when
    a new primer is used, as the padding of the target sequence depends on the
    primer length. It handles both the forward and reverse complement strands
    of the target DNA.

    Attributes:
        padding_len (int): The length of the padding, equal to the primer length.
        primer (Primer): The primer used in this configuration.
        template (DNA): The template DNA sequence.
        template_seq (Dict[DNADirection, str]): A dictionary containing the
                                                padded forward and reverse
                                                template sequences.
        settings (Settings): The settings for replication analysis.
        origin_idx (OriginIdx): The index of replication origins.
    """

    def __init__(
        self,
        template: DNA,
        primer: Primer,
        settings: Settings,
    ) -> None:
        """Construct the Repliconf object.

        Args:
            template (DNA): The template DNA sequence.
            primer (Primer): The primer for this configuration.
            settings (Settings): The replication analysis settings.
        """
        self.padding_len = len(primer)

        self.primer = primer
        self.template = template

        self.template_seq: dict[DNADirection, str] = {}
        # Add padding the 5' end of the template
        self.template_seq[DNADirection.FWD] = template.pad(self.padding_len).seq
        # Add padding to the 3' end of the DNA, compute the complement.
        self.template_seq[DNADirection.REV] = (
            template.reverse().pad(self.padding_len).reverse().complement().seq
        )

        logging.debug(
            f"Repliconf.__init__(): FWD: {self.template_seq[DNADirection.FWD]}"
        )
        logging.debug(
            f"Repliconf.__init__(): REV: {self.template_seq[DNADirection.REV]}"
        )

        self.settings = settings
        self.origin_idx = OriginIdx([], [], False)

    def range(self) -> range:
        """Return the valid search range for replication origins.

        The search range is determined by the length of the padded template
        sequence and the length of the primer.

        Returns:
            range: The valid search range.
        """
        return range(len(self.template_seq[DNADirection.FWD]) - len(self.primer) + 1)

    def slice(self, i: int) -> slice:
        """Return a slice object for constructing a ReplicationOrigin.

        This slice is used to extract the target sequence for a potential
        replication origin at a given index.

        Args:
            i (int): The starting index of the slice.

        Returns:
            slice: A slice object.
        """
        return slice(i, i + len(self.primer))

    def origin(self, direction: DNADirection, i: int) -> ReplicationOrigin:
        """Return the ith ReplicationOrigin.

        This method constructs a `ReplicationOrigin` object for a given
        direction and index.

        Args:
            direction (DNADirection): The direction of the DNA strand.
            i (int): The index of the replication origin.

        Returns:
            ReplicationOrigin: The replication origin object.
        """
        return ReplicationOrigin(
            (
                self.template_seq[direction][self.slice(i)][::-1]
                if direction
                else self.template_seq[direction][self.slice(i)]
            ),
            self.primer.seq[::-1],
            self.settings,
        )

    def search(self) -> None:
        """Search for valid replication origins in both directions."""
        self.origin_idx.clear()
        for direction in [DNADirection.FWD, DNADirection.REV]:
            logging.debug(f"Repliconf.search(): {direction}")
            for i in self.range():
                origin = self.origin(direction, i)
                if (
                    origin.primability > self.settings.primability_cutoff
                    and origin.stability > self.settings.stability_cutoff
                ):
                    logging.debug(
                        f"Repliconf.search(): adding [{direction}, {i}]: {origin}"
                    )
                    if direction:
                        self.origin_idx.fwd.append(i)
                    else:
                        self.origin_idx.rev.append(i)
        self.origin_idx.searched = True

    def __eq__(self, other: object) -> bool:
        """Check if two Repliconf objects are equal.

        Two Repliconf objects are considered equal if they have the same primer
        and template DNA.

        Args:
            other (object): The object to compare with.

        Returns:
            bool: True if the objects are equal, False otherwise.
        """
        if not isinstance(other, Repliconf):
            return NotImplemented
        return self.primer == other.primer and self.template == other.template

    def __hash__(self) -> int:
        """Return a hash value for the Repliconf object.

        The hash is based on the primer and template DNA.

        Returns:
            int: The hash value.
        """
        return hash((self.primer, self.template))

    def __str__(self) -> str:
        """Return the string representation of the Repliconf object.

        Returns:
            str: A string representation of the object.
        """
        return f"ReplicationConfig: Primer: {self.primer}, Target: {self.template}"

    @cached_property
    def amplicon_start(self) -> list[int]:
        """Return the list of amplicon starting positions.

        The starting positions are calculated from the forward replication
        origin indices.

        Returns:
            List[int]: A list of amplicon starting positions.
        """
        return [x - len(self.primer) for x in self.origin_idx.fwd]

    @cached_property
    def amplicon_end(self) -> list[int]:
        """Return the list of amplicon ending positions.

        The ending positions are calculated from the reverse replication
        origin indices.

        Returns:
            List[int]: A list of amplicon ending positions.
        """
        return [x + len(self.primer) for x in self.origin_idx.rev]
