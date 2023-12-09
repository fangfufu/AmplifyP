# -*- coding: utf-8 -*-
"""Amplify P - Replication related."""

from typing import Dict, List, Tuple

from .dna import DNA, Primer, DNADirection
from .origin import ReplicationOrigin
from .settings import Settings

class OriginIndex:
    """The origin index class.

    This class stores the index numbers of the valid replication origins for
    each DNA / DNA direction pair.
    """

    def __init__(self) -> None:
        """Initializes a OriginIndex object."""
        self.__index: Dict[Tuple[DNA, DNADirection], List[int]] = {}

    def __getitem__(self, key: tuple[DNA, DNADirection]) -> List[int]:
        """Return the index number of a valid replication origin."""
        if key not in self.__index:
            return []
        return self.__index[key]

    def append(self, dna: DNA, direction: DNADirection, index: int) -> None:
        """Add a valid index number to the list."""
        if (dna, direction) not in self.__index:
            self.__index[dna, direction] = []
        self.__index[dna, direction].append(index)

    def clear(self, dna: DNA, direction: DNADirection) -> None:
        """Clear the dictionary of the valid index numbers."""
        if (dna, direction) not in self.__index:
            return
        self.__index[dna, direction].clear()

    def clear_all(self) -> None:
        """Clear all the match indices."""
        self.__index.clear()

    def remove(self, dna: DNA, direction: DNADirection, index: int) -> None:
        """Remove the index of a DNA / direction pair."""
        if (dna, direction) not in self.__index:
            return
        self.__index[dna, direction].remove(index)

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

    def __init__(self, template: DNA, primer: Primer, settings: Settings) -> None:
        """Initialize the ReplicationConfig object."""
        min_overlap: int = settings.min_overlap
        padding_len: int = len(primer) - min_overlap

        self.primer = primer
        self.template = template

        # We reverse the sequences here, because Bill's algorithm requires you
        # to count from the 3' end. This is documented in ReplicationOrigin's
        # class docstring.
        self.template_seq: Dict[DNADirection, str] = {}
        self.template_seq[DNADirection.FWD] = template.pad(padding_len).reverse().sequence
        self.template_seq[DNADirection.REV] = (
            template.complement().pad(padding_len).reverse().sequence
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
        return self.primer == other.primer and self.template == other.template

    def __hash__(self) -> int:
        return hash((self.primer, self.template))

    def range(self) -> range:
        """Return the range of the target in ReplicationConfig."""
        return range(
            0, len(self.template_seq[DNADirection.FWD]) - len(self.primer_seq) + 1
        )

    def slice(self, i: int) -> slice:
        """Return the slice of the target in ReplicationConfig."""
        return slice(i, i + len(self.primer_seq))

    def origin(self, direction: DNADirection, i: int) -> ReplicationOrigin:
        """Return the origin of replication."""
        return ReplicationOrigin(
            self.template_seq[direction][self.slice(i)],
            self.primer_seq,
            self.settings,
        )

    # def clear(self) -> None:
    #     """Clear the primer index."""
    #     for direction in DNADirection:
    #         self.primer.index.clear(self.template, direction)

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
                    # self.primer.index.append(
                    #     self.template,
                    #     direction,
                    #     i,
                    # )
                    pass

    def __str__(self) -> str:
        """Return the string representation of the ReplicationConfig object."""
        return f"ReplicationConfig: Primer: {self.primer}, Target: {self.template}"

