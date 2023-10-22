# -*- coding: utf-8 -*-
"""Amplify P - DNA related."""
from enum import Flag, IntEnum, StrEnum


class Nucleotides(StrEnum):
    """Enumeration of valid nucleotide characters for use in AmplifyP."""

    SINGLE = "GATC"
    DOUBLE = "MRWSYK"
    TRIPLE = "VHDB"
    WILDCARD = "N"
    GAP = "-"

    TARGET = SINGLE + WILDCARD
    PRIMER = SINGLE + DOUBLE + TRIPLE + WILDCARD


class DNAType(IntEnum):
    """An enumeration representing the type of DNA."""

    LINEAR = 1
    CIRCULAR = 2
    PRIMER = 3


class DNADirection(Flag):
    """An enumeration representing the direction of DNA."""

    FORWARD = False
    REVERSE = True


class DNA:
    """A class representing a DNA sequence.

    Attributes:
        sequence (str): The DNA sequence.
        dna_type (DNAType): The type of DNA (e.g. primer or target).
        name (str): The name of the DNA object.
        direction (bool | DNADirection): The direction of the DNA sequence.

    Methods:
        lower() -> DNA: Returns the DNA sequence in lower case.
        upper() -> DNA: Returns the DNA sequence in upper case.
        complement() -> DNA: Returns the complement of the DNA sequence.
        reverse() -> DNA: Returns the reverse of the DNA sequence.
        is_complement_of(other: DNA) -> bool: Returns True if the other DNA is
            a complement of this sequence.
        circular_pad() -> DNA: Pads the circular DNA sequence to the required length.
    """

    def __init__(
        self,
        sequence: str,
        dna_type: DNAType = DNAType.LINEAR,
        name: str = " ",
        direction: bool | DNADirection = DNADirection.FORWARD,
    ) -> None:
        """Initializes a DNA object.

        Args:
            sequence (str): The DNA sequence.
            dna_type (DNAType): The type of DNA (e.g. primer or target).
            name (str, optional): The name of the DNA object. Defaults to " ".
            reversed (bool | DNADirection, optional): The direction of the DNA
                sequence. Defaults to DNADirection.FORWARD.

        Raises:
            ValueError: If the DNA sequence contains invalid characters.
        """
        self.__sequence: str = sequence
        self.__dna_type: DNAType = dna_type
        self.__name: str = name
        self.__direction: DNADirection | bool = direction
        check_str = (
            Nucleotides.PRIMER
            if dna_type == DNAType.PRIMER
            else Nucleotides.TARGET + Nucleotides.GAP
        )
        if not set(sequence.upper()) <= set(check_str):
            raise ValueError("DNA sequence contains invalid characters.")

    @property
    def sequence(self) -> str:
        """Return the DNA sequence as a string."""
        return self.__sequence

    @property
    def dna_type(self) -> DNAType:
        """Return the DNA type."""
        return self.__dna_type

    @property
    def name(self) -> str:
        """Return the name of the DNA sequence."""
        return self.__name

    @property
    def direction(self) -> bool | DNADirection:
        """Return the direction of the DNA sequence."""
        return self.__direction

    def lower(self) -> "DNA":
        """Return the DNA sequence in lower case."""
        return DNA(self.sequence.lower(), self.dna_type, self.name, self.direction)

    def upper(self) -> "DNA":
        """Return the DNA sequence in upper case."""
        return DNA(self.sequence.upper(), self.dna_type, self.name, self.direction)

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
            self.direction,
        )

    def reverse(self) -> "DNA":
        """Return the reverse of the DNA sequence."""
        return DNA(self.sequence[::-1], self.dna_type, self.name, not self.direction)

    def __eq__(self, other: object) -> bool:
        """Return True if the DNA sequences are equal."""
        if not isinstance(other, DNA):
            return NotImplemented
        return self.sequence.upper() == other.sequence.upper()

    def is_complement_of(self, other: "DNA") -> bool:
        """Return True if the other DNA is a complement of this sequence."""
        return self.sequence.upper() == other.complement().sequence.upper()

    def __hash__(self) -> int:
        """Return the hash value of the DNA sequence in uppercase."""
        return self.sequence.upper().__hash__()

    def __add__(self, other: object) -> "DNA":
        """Return the concatenation of two DNA sequences."""
        if not isinstance(other, DNA):
            return NotImplemented
        if self.dna_type != other.dna_type:
            return NotImplemented
        return DNA(
            self.sequence + other.sequence, self.dna_type, self.name, self.direction
        )

    def __len__(self) -> int:
        """Return the length of the DNA sequence."""
        return len(self.sequence)

    def circular_pad(self) -> "DNA":
        """Pad the circular DNA sequence to the required length."""
        if self.dna_type == DNAType.CIRCULAR:
            return DNA(
                self.sequence * 3,
                self.dna_type,
                self.name,
                self.direction,
            )
        return NotImplemented

    def __getitem__(self, key: slice) -> "DNA":
        """Return the character at the given index."""
        return DNA(self.sequence[key], self.dna_type, self.name, self.direction)

    def __str__(self) -> str:
        """Return the string representation of the DNA sequence."""
        return self.name

    def __repr__(self) -> str:
        """Return the representation of the DNA sequence."""
        return f"{self.name}: {self.dna_type}, {self.direction}, {self.sequence}"
