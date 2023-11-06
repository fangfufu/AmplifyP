# -*- coding: utf-8 -*-
"""Amplify P - DNA related."""
from enum import Flag, IntEnum, StrEnum
from typing import Dict, List, Tuple


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

    FWD = True
    REV = False

    def __repr__(self) -> str:
        return f"DNADir.{self.name}"

    def __str__(self) -> str:
        return f"{self.name}"


class DNA:
    """A class representing a DNA sequence."""

    def __init__(
        self,
        sequence: str,
        dna_type: DNAType = DNAType.LINEAR,
        name: str | None = None,
        direction: bool | DNADirection = DNADirection.FWD,
    ) -> None:
        """Initializes a DNA object."""
        self.__sequence: str = sequence.strip()
        self.__type: DNAType = dna_type
        if name is None:
            self.__name = sequence
        else:
            self.__name = name.strip()
        self.__direction: DNADirection | bool = direction
        if dna_type == DNAType.CIRCULAR:
            check_str = Nucleotides.CIRCULAR
        elif dna_type == DNAType.LINEAR:
            check_str = Nucleotides.LINEAR
        elif dna_type == DNAType.PRIMER:
            check_str = Nucleotides.PRIMER
        else:
            raise ValueError("Invalid DNA type.")

        if not set(sequence.upper()) <= set(check_str):
            raise ValueError("The DNA sequence contains invalid characters.")

    @property
    def sequence(self) -> str:
        """Return the DNA sequence as a string."""
        return self.__sequence

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
        """
        Return True if the DNA sequences are identical.

        DNA are considered identical if they have the same sequence, direction,
        and type.
        """
        if not isinstance(other, DNA):
            return NotImplemented
        return (
            self.sequence.upper() == other.sequence.upper()
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
        return hash((self.sequence.upper(), self.direction, self.type))

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
        base_str = (
            self.sequence
            if self.direction == DNADirection.FWD
            else self.reverse().sequence
        )
        padding_str = (
            base_str[-i::] if self.type == DNAType.CIRCULAR else Nucleotides.GAP * i
        )
        new_str = padding_str + base_str

        if self.direction == DNADirection.REV:
            new_str = new_str[::-1]

        return DNA(
            new_str,
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


class OriginIndex:
    """The origin index class.

    This class stores the indices of valid replication origins for each DNA /
    DNA direction pair.
    """

    def __init__(self) -> None:
        """Initializes a PrimerIndex object."""
        self.__index: Dict[Tuple[DNA, DNADirection], List[int]] = {}

    def __getitem__(self, key: tuple[DNA, DNADirection]) -> List[int]:
        """Return the match indices of the Primer."""
        if key not in self.__index:
            return []
        return self.__index[key]

    def append(self, dna: DNA, direction: DNADirection, index: int) -> None:
        """Append the match index of the Primer."""
        if (dna, direction) not in self.__index:
            self.__index[dna, direction] = []
        self.__index[dna, direction].append(index)

    def clear(self, dna: DNA, direction: DNADirection) -> None:
        """Clear the match index of the Primer."""
        if (dna, direction) not in self.__index:
            return
        self.__index[dna, direction].clear()

    def remove(self, dna: DNA, direction: DNADirection, index: int) -> None:
        """Remove the match index of the Primer."""
        if (dna, direction) not in self.__index:
            return
        self.__index[dna, direction].remove(index)


class Primer(DNA):
    """
    A class representing a Primer sequence.

    Primer is a subclass of DNA, but it has an additional attribute called
    index, which is a dictionary that stores the valid replication origin for
    each DNA sequence / DNA direction pair.
    """

    def __init__(
        self,
        sequence: str,
        name: str | None = None,
    ) -> None:
        """Initializes a Primer object."""
        super().__init__(sequence, DNAType.PRIMER, name, DNADirection.FWD)
        self.__index: OriginIndex = OriginIndex()

    @property
    def index(self) -> OriginIndex:
        """Return the match indices of the Primer."""
        return self.__index
