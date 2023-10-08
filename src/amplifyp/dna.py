# -*- coding: utf-8 -*-
"""Amplify P - DNA related."""
from dataclasses import dataclass
from enum import IntEnum, StrEnum


class Nucleotides(StrEnum):
    """Enumeration of valid nucleotide characters for use in AmplifyP."""

    SINGLE = "GATC"
    DOUBLE = "MRWSYK"
    TRIPLE = "VHDB"
    WILDCARD = "N"

    TARGET = SINGLE + WILDCARD
    PRIMER = TARGET + DOUBLE + TRIPLE


class DNAType(IntEnum):
    """An enumeration representing the type of DNA."""

    DEFAULT = 0
    PRIMER = 1
    CIRCULAR = 2


@dataclass(slots=True)
class DNA:
    """A class representing a DNA sequence.

    Attributes:
        sequence (str): The DNA sequence.
        name (str): The name of the DNA sequence (optional).
        dna_type (DNAType): The type of the DNA sequence (default is DNAType.DEFAULT).
    """

    sequence: str
    name: str = ""
    dna_type: DNAType = DNAType.DEFAULT

    def __post_init__(self) -> None:
        """Validate the DNA sequence.

        Raises:
            ValueError: If the DNA sequence contains invalid characters.
        """
        check_str = (
            Nucleotides.PRIMER
            if self.dna_type == DNAType.PRIMER
            else Nucleotides.TARGET
        )
        if not set(self.sequence.upper()) <= set(check_str):
            raise ValueError("DNA sequence contains invalid characters.")

    def lower(self) -> "DNA":
        """Return the DNA sequence in lower case."""
        return DNA(self.sequence.lower(), self.name, self.dna_type)

    def upper(self) -> "DNA":
        """Return the DNA sequence in upper case."""
        return DNA(self.sequence.upper(), self.name, self.dna_type)

    def complement(self) -> "DNA":
        """Return the complement of the DNA sequence.

        Note that the complement of non-ACGT bases are undefined.
        """
        return DNA(
            self.sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1],
            self.name,
            self.dna_type,
        )

    def __eq__(self, other: object) -> bool:
        """Return True if the DNA sequences are equal."""
        if not isinstance(other, DNA):
            return NotImplemented

        self_seq = self.upper().sequence
        other_seq = other.upper().sequence
        other_complement_seq = other.upper().complement().sequence
        return self_seq in [other_seq, other_complement_seq]

    def __add__(self, other: object) -> "DNA":
        """Return the concatenation of two DNA sequences."""
        if not isinstance(other, DNA):
            return NotImplemented

        return DNA(self.sequence + other.sequence, self.name, self.dna_type)
