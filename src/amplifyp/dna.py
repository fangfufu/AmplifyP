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
    BLANK = "-"

    TARGET = SINGLE + WILDCARD
    PRIMER = TARGET + DOUBLE + TRIPLE


class DNAType(IntEnum):
    """An enumeration representing the type of DNA."""

    LINEAR = 1
    CIRCULAR = 2
    PRIMER = 3


@dataclass(slots=True)
class DNA:
    """A class representing a DNA sequence.

    Attributes:
        sequence (str): The DNA sequence.
        dna_type (DNAType): The type of the DNA sequence (default is DNAType.LINEAR).
        name (str): The name of the DNA sequence (optional).
    """

    sequence: str
    dna_type: DNAType = DNAType.LINEAR
    name: str = ""
    reversed: bool = False

    def __post_init__(self) -> None:
        """Validate the DNA sequence.

        Raises:
            ValueError: If the DNA sequence contains invalid characters.
        """
        check_str = (
            Nucleotides.PRIMER
            if self.dna_type == DNAType.PRIMER
            else Nucleotides.TARGET + Nucleotides.BLANK
        )
        if not set(self.sequence.upper()) <= set(check_str):
            raise ValueError("DNA sequence contains invalid characters.")

    def lower(self) -> "DNA":
        """Return the DNA sequence in lower case."""
        return DNA(self.sequence.lower(), self.dna_type, self.name, self.reversed)

    def upper(self) -> "DNA":
        """Return the DNA sequence in upper case."""
        return DNA(self.sequence.upper(), self.dna_type, self.name, self.reversed)

    def complement(self) -> "DNA":
        """Return the complement of the DNA sequence.

        Please refer to the following link:
        https://en.wikipedia.org/wiki/Nucleic_acid_notation
        """
        return DNA(
            self.sequence.translate(
                str.maketrans("ACGTMKRYBDHVacgtmkrybdhv", "TGCAKMYRVHDBtgcakmyrvhdb")
            )[::-1],
            self.dna_type,
            self.name,
            self.reversed,
        )

    def reverse(self) -> "DNA":
        """Return the reverse of the DNA sequence."""
        return DNA(self.sequence[::-1], self.dna_type, self.name, not self.reversed)

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
        if self.dna_type != other.dna_type:
            return NotImplemented
        return DNA(
            self.sequence + other.sequence, self.dna_type, self.name, self.reversed
        )

    def __len__(self) -> int:
        """Return the length of the DNA sequence."""
        return len(self.sequence)

    def pad(self, size: int) -> "DNA":
        """Pad the DNA sequence with blank characters."""
        if self.dna_type == DNAType.LINEAR:
            return DNA(
                Nucleotides.BLANK * size + self.sequence + Nucleotides.BLANK * size,
                self.dna_type,
                self.name,
                self.reversed,
            )
        if self.dna_type == DNAType.CIRCULAR:
            return DNA(
                self.sequence[-size:] + self.sequence + self.sequence[:size],
                self.dna_type,
                self.name,
                self.reversed,
            )
        return NotImplemented

    def __getitem__(self, key: slice) -> "DNA":
        """Return the character at the given index."""
        return DNA(self.sequence[key], self.dna_type, self.name, self.reversed)
