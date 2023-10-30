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

    TARGET = SINGLE + WILDCARD + GAP
    PRIMER = SINGLE + DOUBLE + TRIPLE + WILDCARD + GAP


class DNAType(IntEnum):
    """An enumeration representing the type of DNA."""

    LINEAR = 1
    CIRCULAR = 2
    PRIMER = 3

    def __repr__(self) -> str:
        return f"DNAType.{self.name}"

    def __str__(self) -> str:
        return f"{self.name}"


class DNADirection(Flag):
    """
    An enumeration representing the direction of DNA.

    Forward DNA (True) means DNA going from 5' to 3'.
    Reverse DNA (False) means DNA going from 3' to 5'.
    """

    FORWARD = True
    REVERSE = False

    def __repr__(self) -> str:
        return f"DNADirection.{self.name}"

    def __str__(self) -> str:
        return f"{self.name}"


class DNA:
    """A class representing a DNA sequence."""

    def __init__(
        self,
        sequence: str,
        dna_type: DNAType = DNAType.LINEAR,
        name: str | None = None,
        direction: bool | DNADirection = DNADirection.FORWARD,
    ) -> None:
        """Initializes a DNA object."""
        self.__sequence: str = sequence.strip()
        self.__type: DNAType = dna_type
        if name is None:
            self.__name = sequence
        else:
            self.__name = name.strip()
        self.__direction: DNADirection | bool = direction
        check_str = (
            Nucleotides.PRIMER if dna_type == DNAType.PRIMER else Nucleotides.TARGET
        )
        if not set(sequence.upper()) <= set(check_str):
            raise ValueError("DNA sequence contains invalid characters.")

    @property
    def sequence(self) -> str:
        """Return the DNA sequence as a string."""
        return self.__sequence

    @sequence.setter
    def sequence(self, value: str) -> None:
        """Set the DNA sequence."""
        self.__sequence = value.strip()

    @property
    def type(self) -> DNAType:
        """Return the DNA type."""
        return self.__type

    @type.setter
    def type(self, value: DNAType) -> None:
        """Set the DNA type."""
        self.__type = value

    @property
    def name(self) -> str:
        """Return the name of the DNA sequence."""
        return self.__name

    @name.setter
    def name(self, value: str) -> None:
        """Set the name of the DNA sequence."""
        self.__name = value.strip()

    @property
    def direction(self) -> bool | DNADirection:
        """Return the direction of the DNA sequence."""
        return self.__direction

    @direction.setter
    def direction(self, value: bool | DNADirection) -> None:
        """Set the direction of the DNA sequence."""
        self.__direction = value

    def lower(self) -> "DNA":
        """Return the DNA sequence in lower case."""
        return DNA(self.sequence.lower(), self.type, self.name, self.direction)

    def upper(self) -> "DNA":
        """Return the DNA sequence in upper case."""
        return DNA(self.sequence.upper(), self.type, self.name, self.direction)

    def complement(self) -> "DNA":
        """Return the complement of the DNA sequence."""
        return DNA(
            self.sequence.translate(
                str.maketrans("ACGTMKRYBDHVacgtmkrybdhv", "TGCAKMYRVHDBtgcakmyrvhdb")
            ),
            self.type,
            self.name,
            not self.direction,
        )

    def reverse(self) -> "DNA":
        """Return the reverse of the DNA sequence."""
        return DNA(self.sequence[::-1], self.type, self.name, not self.direction)

    def __eq__(self, other: object) -> bool:
        """Return True if the DNA sequences are identical."""
        if not isinstance(other, DNA):
            return NotImplemented
        return (
            self.sequence.upper() == other.sequence.upper()
            and self.direction == other.direction
        )

    def is_complement_of(self, other: "DNA") -> bool:
        """Return True if the other DNA is a complement of this sequence."""
        return (
            self.sequence.upper() == other.complement().sequence.upper()
            and self.direction != other.direction
        )

    def __len__(self) -> int:
        """Return the length of the DNA sequence."""
        return len(self.sequence)

    def pad(self, i: int) -> "DNA":
        """Pad the DNA sequence to the required length."""
        if self.direction == DNADirection.REVERSE:
            new_base_str = self.reverse().sequence
        else:
            new_base_str = self.sequence

        if self.type == DNAType.CIRCULAR:
            padding_str = new_base_str[-i::]
        else:
            padding_str = Nucleotides.GAP * i

        new_return_str = padding_str + new_base_str

        if self.direction == DNADirection.REVERSE:
            new_return_str = new_return_str[::-1]

        return DNA(
            new_return_str,
            self.type,
            self.name,
            self.direction,
        )

    def __getitem__(self, key: slice) -> "DNA":
        """Return the nucleotides at the given index."""
        return DNA(self.sequence[key], self.type, self.name, self.direction)

    def __str__(self) -> str:
        """Return the string representation of the DNA sequence."""
        return f"DNA: {self.name}"

    def __repr__(self) -> str:
        """Return the representation of the DNA sequence."""
        return (
            f"DNA: Name: {self.name}, {repr(self.type)}, "
            + f"{repr(self.direction)}, Seq: {self.sequence}"
        )
