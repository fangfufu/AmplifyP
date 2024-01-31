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

    CIRCULAR = SINGLE + WILDCARD
    LINEAR = CIRCULAR + GAP
    PRIMER = SINGLE + DOUBLE + TRIPLE + WILDCARD


class DNAType(IntEnum):
    """An enumeration representing the type of DNA."""

    LINEAR = 1
    CIRCULAR = 2
    PRIMER = 3


class DNADirection(Flag):
    """
    An enumeration representing the direction of DNA.

    Forward DNA (True) means DNA going from 5' to 3'.
    Reverse DNA (False) means DNA going from 3' to 5'.
    """

    FWD = True
    REV = False


class DNA:
    """A class representing a DNA sequence."""

    def __init__(
        self,
        seq: str,
        dna_type: DNAType = DNAType.LINEAR,
        name: str | None = None,
        direction: bool | DNADirection = DNADirection.FWD,
    ) -> None:
        """Initializes a DNA object."""
        self.__seq: str = "".join(seq.split())
        self.__type: DNAType = dna_type

        if name is None:
            self.__name = seq
        else:
            self.__name = name.strip()
        self.__direction: bool | DNADirection = direction

        if dna_type == DNAType.LINEAR:
            check_str = Nucleotides.LINEAR
        elif dna_type == DNAType.CIRCULAR:
            check_str = Nucleotides.CIRCULAR
        elif dna_type == DNAType.PRIMER:
            check_str = Nucleotides.PRIMER
        else:
            raise TypeError("Invalid DNA type.")

        if not set(self.__seq.upper()) <= set(check_str):
            raise ValueError("The DNA sequence contains invalid characters.")

    @property
    def seq(self) -> str:
        """Return the DNA sequence as a string."""
        return self.__seq

    @property
    def type(self) -> DNAType:
        """Return the DNA type."""
        return self.__type

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

    def lower(self) -> "DNA":
        """Return the DNA sequence in lower case."""
        return DNA(self.seq.lower(), self.type, self.name, self.direction)

    def upper(self) -> "DNA":
        """Return the DNA sequence in upper case."""
        return DNA(self.seq.upper(), self.type, self.name, self.direction)

    def complement(self) -> "DNA":
        """Return the complement of the DNA sequence."""
        return DNA(
            self.seq.translate(
                str.maketrans("ACGTMKRYBDHVacgtmkrybdhv", "TGCAKMYRVHDBtgcakmyrvhdb")
            ),
            self.type,
            self.name,
            not self.direction,
        )

    def reverse(self) -> "DNA":
        """Return the reverse of the DNA sequence."""
        return DNA(self.seq[::-1], self.type, self.name, not self.direction)

    def __eq__(self, other: object) -> bool:
        """
        Return True if the DNA sequences are identical.

        DNA are considered identical if they have the same sequence, direction,
        and type.
        """
        if not isinstance(other, DNA):
            return NotImplemented
        return (
            self.seq.upper() == other.seq.upper()
            and self.direction == other.direction
            and self.type == other.type
        )

    def __hash__(self) -> int:
        """
        Return a hash value for the DNA object.

        The hash value is calculated based on the DNA sequence, direction,
        and type.

        Returns:
            int: The hash value for the DNA object.
        """
        return hash((self.seq.upper(), self.direction, self.type))

    def is_complement_of(self, other: "DNA") -> bool:
        """Return True if the other DNA is a complement of this sequence."""
        return (
            self.seq.upper() == other.complement().seq.upper()
            and self.direction != other.direction
        )

    def __len__(self) -> int:
        """Return the length of the DNA sequence."""
        return len(self.seq)

    def pad(self, i: int) -> "DNA":
        """Add padding to beginning of the DNA sequence string."""

        if self.type == DNAType.LINEAR:
            padding_str = Nucleotides.GAP * i
        elif self.type == DNAType.CIRCULAR:
            padding_str = self.seq[-i::]
        else:
            raise TypeError("Invalid DNA type for padding operation.")

        new_str = padding_str + self.seq

        return DNA(
            new_str,
            self.type,
            self.name,
            self.direction,
        )

    def rot(self, i: int) -> "DNA":
        """Rotate the DNA sequence by i bases."""
        if self.type == DNAType.LINEAR:
            raise TypeError("Rotation is unsupported for linear DNA.")
        padded_dna = self.pad(i)

        return padded_dna[0 : len(self.seq)]

    def __getitem__(self, key: slice) -> "DNA":
        """Return the nucleotides at the given index."""
        return DNA(self.seq[key], self.type, self.name, self.direction)

    def __str__(self) -> str:
        """Return the representation of the DNA sequence."""
        return (
            f"DNA: {self.name}, {self.type.name}, {DNADirection(self.direction).name}"
        )


class Primer(DNA):
    """
    A class representing a Primer sequence.

    Primer is a subclass of DNA
    """

    def __init__(
        self,
        sequence: str,
        name: str | None = None,
    ) -> None:
        """Initializes a Primer object."""
        super().__init__(sequence, DNAType.PRIMER, name, DNADirection.FWD)
