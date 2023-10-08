# -*- coding: utf-8 -*-
"""Amplify P - Data types."""
from dataclasses import dataclass
from enum import IntEnum, StrEnum
from typing import Iterator, List, Tuple


class DNAType(IntEnum):
    """An enumeration representing the type of DNA."""

    DEFAULT = 0
    PRIMER = 1
    PLASMID = 2


class Nucleotides(StrEnum):
    """Enumeration of valid nucleotide characters for use in AmplifyP."""

    SINGLE = "GATC"
    DOUBLE = "MRWSYK"
    TRIPLE = "VHDB"
    WILDCARD = "N"

    TARGET = SINGLE + WILDCARD
    PRIMER = TARGET + DOUBLE + TRIPLE


class DNA:
    """A class representing a DNA sequence."""

    def __init__(
        self,
        sequence: str,
        dna_type: DNAType = DNAType.DEFAULT,
        name: str = "",
    ) -> None:
        """Construct a DNA sequence.

        Args:
            sequence (str): The DNA sequence to be constructed.
            primer (bool, optional): Whether the sequence is a primer or not.
                Defaults to False.
            name (str, optional): The name of the sequence. Defaults to "".

        Raises:
            ValueError: If the DNA sequence contains invalid characters.
        """
        check_str = (
            Nucleotides.PRIMER if dna_type == DNAType.PRIMER else Nucleotides.TARGET
        )
        if not set(sequence.upper()) <= set(check_str):
            raise ValueError("DNA sequence contains invalid characters.")
        self._sequence = sequence
        self._type = dna_type
        self._name = name

    def __str__(self) -> str:
        """Return the string representation of a DNA sequence."""
        return self.sequence

    def __repr__(self) -> str:
        """Return the string representation of a DNA sequence."""
        return self.sequence

    def __eq__(self, other: object) -> bool:
        """Return True if the DNA sequences are equal."""
        if isinstance(other, DNA):
            return self.sequence.upper() == other.sequence.upper()
        return NotImplemented

    @property
    def name(self) -> str:
        """Return the name of the DNA sequence."""
        return self._name

    @name.setter
    def name(self, value: str) -> None:
        """Set the name of the DNA sequence."""
        self._name = value

    @property
    def type(self) -> DNAType:
        """Return True if the DNA sequence is a primer."""
        return self._type

    @property
    def sequence(self) -> str:
        """Return the DNA sequence."""
        return self._sequence

    def reverse(self) -> str:
        """Return the reverse complement of the DNA sequence."""
        return self._sequence[::-1]

    def lower(self) -> str:
        """Return the DNA sequence in lower case."""
        return self._sequence.lower()

    def upper(self) -> str:
        """Return the DNA sequence in upper case."""
        return self._sequence.upper()

    def complement(self) -> str:
        """Return the complement of the DNA sequence.

        Note that the complement of non-ACGT bases are undefined.
        """
        return self._sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))


class LengthWiseWeightTbl:
    """A class representing a run-length weight table."""

    def __init__(
        self,
        size: int,
        init_weight: float = 0,
        overrides: List[Tuple[int, float]] | None = None,
    ) -> None:
        """Initialize a Run-length Weight Table.

        Args:
            size (int): The size of the weight table.
            init_weight (float, optional): The initial weight for all elements
                in the table. Defaults to 0.
            overrides (List[Tuple[int, float]] | None, optional): A list
                of (index, weight) tuples to override the initial weights.
                Defaults to None.
        """
        self._weight = [init_weight] * size

        if overrides is not None:
            for item in overrides:
                key, value = item
                self._weight[key] = value

    def __getitem__(self, key: int) -> float:
        """Return the weight of at certain run-length."""
        return self._weight[key]

    def __setitem__(self, key: int, value: float) -> None:
        """Set the weight at a certain run-length."""
        self._weight[key] = value

    def __len__(self) -> int:
        """Return the size of the Run-length Weight Table."""
        return len(self._weight)

    def __str__(self) -> str:
        """Return the string representation of the Run-length Weight Table."""
        return str(self._weight)

    def __repr__(self) -> str:
        """Return the string representation of the Run-length Weight Table."""
        return str(self._weight)

    def __iter__(self) -> Iterator[float]:
        """Return the iterator of the Run-length Weight Table."""
        return iter(self._weight)


class BasePairWeights:
    """Nucleotide Pairwise Weight Table."""

    def __init__(self, row: str, column: str, weight: List[List[float]]) -> None:
        """Construct a Nucleotide Pairwise Weight Table.

        Args:
            row (str): A string representing the row labels of the weight table.
            column (str): A string representing the column labels of the weight
                table.
            weight (List[List[float]]): A 2D list of floats representing the
                weights of each nucleotide pair.

        Raises:
            ValueError: If the length of the weight table does not match the
                length of the row or column labels.
        """
        self._row = dict(zip(list(row), range(len(row))))
        self._column = dict(zip(list(column), range(len(column))))
        if len(weight) != len(row):
            raise ValueError(
                "NucleotidePairwiseWeightTbl: row length mismatch at initialisation."
            )
        for i in range(len(row)):
            if len(weight[i]) != len(column):
                raise ValueError(
                    "NucleotidePairwiseWeightTbl: column length mismatch at initialisation."
                )
        self._weight = weight

    @property
    def row(self) -> list[str]:
        """Return the row nucleotides."""
        return list(self._row.keys())

    @property
    def column(self) -> list[str]:
        """Return the column nucleotides."""
        return list(self._column.keys())

    def __getitem__(self, key: tuple[str, str]) -> float:
        """Return the weight of at certain nucleotide pair."""
        row, column = key
        return self._weight[self._row[row]][self._column[column]]

    def __setitem__(self, key: tuple[str, str], value: float) -> None:
        """Set the weight at a certain nucleotide pair."""
        row, column = key
        self._weight[self._row[row]][self._column[column]] = value

    def __len__(self) -> int:
        """Return the size of the Run-length Weight table."""
        return len(self.row) * len(self.column)

    def __str__(self) -> str:
        """Return the string representation of the table."""
        return str(self._weight)

    def __repr__(self) -> str:
        """Return the string representation of the table."""
        return str(self._weight)


@dataclass(frozen=True)
class Replisome:
    """A class representing a replisome.

    Attributes:
        target (DNA): The DNA sequence being replicated.
        primer (DNA): The DNA sequence that serves as a starting point for DNA
            synthesis.
        base_pair_scores (BasePairWeights): The weights assigned to each base
            pair in the DNA sequence.
        match_weight (LengthWiseWeightTbl): The weights assigned to matches
            between the target and primer sequences.
        run_weight (LengthWiseWeightTbl): The weights assigned to runs of
            consecutive matches between the target and primer sequences.
    """

    target: DNA
    primer: DNA
    base_pair_scores: BasePairWeights
    match_weight: LengthWiseWeightTbl
    run_weight: LengthWiseWeightTbl
