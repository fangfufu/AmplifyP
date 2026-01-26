# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""DNA-related classes and enumerations for AmplifyP."""

from __future__ import annotations

from enum import Flag, IntEnum, StrEnum

GLOBAL_COMPLEMENT_TABLE = str.maketrans(
    "ACGTMKRYBDHVacgtmkrybdhv", "TGCAKMYRVHDBtgcakmyrvhdb"
)


class Nucleotides(StrEnum):
    """Enumeration of valid nucleotide characters for use in AmplifyP.

    This class defines groups of nucleotide characters that are valid in
    different contexts within the AmplifyP application. It includes single
    nucleotides, ambiguous (double and triple) nucleotides, wildcards, and gaps.
    It also defines composite groups for circular, linear, and primer DNA types.

    Attributes:
        SINGLE: Standard single nucleotides (G, A, T, C).
        DOUBLE: Ambiguous nucleotides representing two possibilities (M, R, W,
            S, Y, K).
        TRIPLE: Ambiguous nucleotides representing three possibilities (V, H, D,
            B).
        WILDCARD: Wildcard character representing any nucleotide (N).
        GAP: Gap character (-).
        PRIMER: Valid characters for primers (SINGLE + DOUBLE + TRIPLE +
            WILDCARD).
    """

    SINGLE = "GATC"
    DOUBLE = "MRWSYK"
    TRIPLE = "VHDB"
    WILDCARD = "N"
    GAP = "-"
    TEMPLATE = SINGLE + WILDCARD + GAP
    PRIMER = SINGLE + DOUBLE + TRIPLE + WILDCARD


class DNAType(IntEnum):
    """An enumeration representing the type of DNA.

    This enumeration is used to specify whether a DNA sequence is linear,
    circular, or a primer. This distinction is important for certain
    operations, such as padding and rotation.

    Attributes:
        LINEAR: Represents linear DNA (e.g., chromosomes).
        CIRCULAR: Represents circular DNA (e.g., plasmids).
        PRIMER: Represents a primer sequence.
    """

    LINEAR = 1
    CIRCULAR = 2
    PRIMER = 3


class DNADirection(Flag):
    """An enumeration representing the direction of a DNA strand.

    This enumeration defines the two possible directions for a DNA strand:
    forward (5' to 3') and reverse (3' to 5'). It is used to track the
    orientation of DNA sequences.

    Attributes:
        FWD: Represents the forward direction (5' -> 3').
        REV: Represents the reverse direction (3' -> 5').
    """

    FWD = True
    REV = False

    def __str__(self) -> str:
        """The direction as a string."""
        return "Forward" if self == DNADirection.FWD else "Reverse"


class DNA:
    """A class representing a DNA sequence.

    This class encapsulates a DNA sequence and its properties, such as its type
    (linear, circular, or primer), name, and direction. It provides methods for
    manipulating and analyzing the DNA sequence, such as complementing,
    reversing, and padding.

    Attributes:
        seq (str): The DNA sequence string.
        type (DNAType): The type of the DNA sequence (LINEAR, CIRCULAR, or
            PRIMER).
        name (str): The name identifier of the DNA sequence.
        direction (bool | DNADirection): The direction of the DNA sequence (FWD
            or REV).
    """

    def __init__(
        self,
        seq: str,
        dna_type: DNAType = DNAType.LINEAR,
        name: str | None = None,
        direction: bool | DNADirection = DNADirection.FWD,
    ) -> None:
        """Initializes a new DNA object.

        Args:
            seq (str): The raw DNA sequence. Whitespace will be removed.
            dna_type (DNAType, optional): The type of DNA (LINEAR, CIRCULAR,
                PRIMER). Defaults to DNAType.LINEAR.
            name (str, optional): The name of the DNA sequence. If None, the
                sequence itself is used as the name. Defaults to None.
            direction (bool | DNADirection, optional): The direction of the DNA
                sequence. Defaults to DNADirection.FWD.

        Raises:
            TypeError: If the provided `dna_type` is invalid.
            ValueError: If the `seq` contains characters invalid for the
                specified `dna_type`.
        """
        self.__seq: str = "".join(seq.split())
        self.__type: DNAType = dna_type

        if name is None:
            self.__name = seq
        else:
            self.__name = name.strip()
        self.__direction: bool | DNADirection = direction

        if dna_type == DNAType.LINEAR or dna_type == DNAType.CIRCULAR:
            check_str = Nucleotides.TEMPLATE
        elif dna_type == DNAType.PRIMER:
            check_str = Nucleotides.PRIMER
        else:
            raise TypeError("Invalid DNA type.")

        # Optimization: cache the uppercase sequence to avoid repeated
        # allocations in _count_bases and other methods.
        seq_upper = self.__seq.upper()
        if seq_upper == self.__seq:
            # If the sequence is already uppercase, avoid duplicating the
            # string.
            self._seq_upper = self.__seq
        else:
            self._seq_upper = seq_upper

        invalid_chars = set(self._seq_upper) - set(check_str)
        if invalid_chars:
            raise ValueError(
                f"The DNA sequence contains invalid characters: {
                    ', '.join(sorted(invalid_chars))
                }"
            )

    @property
    def seq(self) -> str:
        """The DNA sequence as a string."""
        return self.__seq

    @property
    def type(self) -> DNAType:
        """The DNA type (LINEAR, CIRCULAR, or PRIMER)."""
        return self.__type

    @property
    def name(self) -> str:
        """The name of the DNA sequence."""
        return self.__name

    @name.setter
    def name(self, value: str) -> None:
        """Set the name of the DNA sequence.

        Args:
            value (str): The new name. Whitespace will be stripped.
        """
        self.__name = value.strip()

    @property
    def direction(self) -> bool | DNADirection:
        """The direction of the DNA sequence."""
        return self.__direction

    def lower(self) -> DNA:
        """Return a copy of the DNA object with the sequence in lowercase.

        Returns:
            DNA: A new DNA object with the lowercase sequence.
        """
        return DNA(self.seq.lower(), self.type, self.name, self.direction)

    def upper(self) -> DNA:
        """Return a copy of the DNA object with the sequence in uppercase.

        Returns:
            DNA: A new DNA object with the uppercase sequence.
        """
        return DNA(self.seq.upper(), self.type, self.name, self.direction)

    def complement(self) -> DNA:
        """Return the complement of the DNA sequence.

        The complement is calculated by swapping A<->T, C<->G, etc.
        The direction of the returned DNA object is inverted relative to the
        original.

        Returns:
            DNA: A new DNA object representing the complement sequence.
        """
        return DNA(
            self.seq.translate(GLOBAL_COMPLEMENT_TABLE),
            self.type,
            self.name,
            not self.direction,
        )

    def reverse(self) -> DNA:
        """Return the reverse of the DNA sequence.

        The direction of the returned DNA object is inverted relative to the
        original.

        Returns:
            DNA: A new DNA object representing the reversed sequence.
        """
        return DNA(self.seq[::-1], self.type, self.name, not self.direction)

    def reverse_complement(self) -> DNA:
        """Return the reverse complement of the DNA sequence.

        This is equivalent to calling `.reverse().complement()`.
        The direction of the returned DNA object is inverted relative to the
        original.

        Returns:
            DNA: A new DNA object representing the reverse complement sequence.
        """
        return self.reverse().complement()

    def __eq__(self, other: object) -> bool:
        """Check if two DNA objects are equal.

        Equality is determined by case-insensitive sequence comparison,
        matching direction, and matching DNA type.

        Args:
            other (object): The object to compare with.

        Returns:
            bool: True if the objects are equal, False otherwise.
        """
        if not isinstance(other, DNA):
            return NotImplemented
        return (
            self._seq_upper == other._seq_upper
            and self.direction == other.direction
            and self.type == other.type
        )

    def __hash__(self) -> int:
        """Return a hash value for the DNA object.

        The hash is computed from the uppercase sequence, direction, and type.

        Returns:
            int: The hash value.
        """
        return hash((self._seq_upper, self.direction, self.type))

    def __len__(self) -> int:
        """Return the length of the DNA sequence.

        Returns:
            int: The length of the sequence.
        """
        return len(self.seq)

    def __getitem__(self, key: slice) -> DNA:
        """Get a slice of the DNA sequence.

        Args:
            key (slice): The slice object specifying the range.

        Returns:
            DNA: A new DNA object representing the subsequence.
        """
        return DNA(self.seq[key], self.type, self.name, self.direction)

    def __str__(self) -> str:
        """Return a string representation of the DNA object.

        Returns:
            str: A formatted string containing name, type, and direction.
        """
        return (
            f"DNA: {self.name}, {self.type.name}, "
            f"{DNADirection(self.direction).name}"
        )

    def __add__(self, other: DNA) -> DNA:
        """Concatenate two DNA sequences.

        Args:
            other (DNA): The DNA sequence to append.

        Returns:
            DNA: A new DNA object with the combined sequence. The type is set
                to LINEAR.
        """
        return DNA(
            self.seq + other.seq, DNAType.LINEAR, self.name, self.direction
        )

    def _count_bases(self, bases: str) -> int:
        """Count the occurrences of any character from `bases` in the sequence.

        Helper method to sum the counts of individual characters provided in the
        `bases` argument within the DNA sequence.

        Args:
            bases (str): A string containing the characters to be counted.

        Returns:
            int: The total count of the specified characters in the sequence.
        """
        return sum(self._seq_upper.count(base) for base in bases)

    def count_at(self) -> int:
        """Count the number of A, T, or W (A/T ambiguous) bases.

        Delegates to `_count_bases` to count "A", "T", and "W".

        Returns:
            int: The total count of A, T, and W bases in the sequence.
        """
        return self._count_bases("ATW")

    def count_cg(self) -> int:
        """Count the number of C, G, or S (C/G ambiguous) bases.

        Delegates to `_count_bases` to count "C", "G", and "S".

        Returns:
            int: The total count of C, G, and S bases in the sequence.
        """
        return self._count_bases("CGS")

    def ratio_at(self) -> float:
        """Calculate the ratio of A, T, or W bases in the sequence.

        Returns:
            float: The ratio of A/T/W bases. Returns 0.0 if sequence is empty.
        """
        if len(self) == 0:
            return 0.0
        return self.count_at() / len(self)

    def ratio_cg(self) -> float:
        """Calculate the ratio of C, G, or S bases in the sequence.

        Returns:
            float: The ratio of C/G/S bases. Returns 0.0 if sequence is empty.
        """
        if len(self) == 0:
            return 0.0
        return self.count_cg() / len(self)


class Primer(DNA):
    """A class representing a primer sequence.

    A primer is a short single-stranded DNA sequence used as a starting point
    for DNA synthesis. This class inherits from `DNA` but enforces
    `DNAType.PRIMER` and `DNADirection.FWD`.
    """

    def __init__(
        self,
        sequence: str,
        name: str | None = None,
    ) -> None:
        """Initializes a Primer object.

        Args:
            sequence (str): The nucleotide sequence of the primer.
            name (str, optional): An optional name for the primer. Defaults to
                None.
        """
        super().__init__(sequence, DNAType.PRIMER, name, DNADirection.FWD)
