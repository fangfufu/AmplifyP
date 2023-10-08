# -*- coding: utf-8 -*-
"""Amplify P - settings related.

Don't ask me how these values were chosen, I copied them from
https://github.com/wrengels/Amplify4/blob/master/Amplify4/Amplify4/TargDelegate.swift
"""

from typing import Final, Iterator, List, Tuple
from .dna import Nucleotides


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


class BasePairWeightsTbl:
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
        if Nucleotides.BLANK in key:
            return 0
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


DEFAULT_MATCH_WEIGHTS: Final[LengthWiseWeightTbl] = LengthWiseWeightTbl(
    size=200,
    init_weight=1,
    overrides=[
        (0, 30),
        (1, 20),
        (2, 10),
        (3, 10),
        (4, 9),
        (5, 9),
        (6, 8),
        (7, 7),
        (8, 6),
        (9, 5),
        (10, 5),
        (11, 4),
        (12, 3),
        (13, 2),
        (14, 1),
    ],
)
DEFAULT_RUN_WEIGHTS: Final[LengthWiseWeightTbl] = LengthWiseWeightTbl(
    size=200,
    init_weight=186,
    overrides=[(0, 100), (1, 150), (2, 175), (3, 182), (4, 186)],
)

DEFAULT_BASE_PAIR_WEIGHTS: Final[BasePairWeightsTbl] = BasePairWeightsTbl(
    row=Nucleotides.PRIMER,
    column=Nucleotides.TARGET,
    weight=[
        [100, 0, 0, 0, 30],
        [0, 100, 0, 0, 30],
        [0, 0, 100, 0, 30],
        [0, 0, 0, 100, 30],
        [0, 70, 0, 70, 30],
        [70, 70, 0, 0, 30],
        [0, 70, 70, 0, 30],
        [70, 0, 0, 70, 30],
        [0, 0, 70, 70, 30],
        [70, 0, 70, 0, 30],
        [50, 50, 0, 50, 30],
        [0, 50, 50, 50, 30],
        [50, 50, 50, 0, 30],
        [50, 0, 50, 50, 30],
        [30, 30, 30, 30, 30],
    ],
)
DEFAULT_PRIMER_DIMER_WEIGHTS: Final[BasePairWeightsTbl] = BasePairWeightsTbl(
    row=Nucleotides.PRIMER,
    column=Nucleotides.PRIMER,
    weight=[
        [-20, -20, -20, 30, 5, -20, -20, 5, 5, -20, -3, -3, -20, -3, -8],
        [-20, -20, 20, -20, -20, -20, 0, -20, 0, 0, -20, -7, -7, -7, -10],
        [-20, 20, -20, -20, 0, 0, 0, -20, -20, -20, -7, -7, -7, -20, -10],
        [30, -20, -20, -20, -20, 5, -20, 5, -20, 5, -3, -20, -3, -3, -8],
        [5, -20, 0, -20, -20, -8, -10, -8, -10, 3, -12, -13, -5, -5, -9],
        [-20, -20, 0, 5, -8, -20, -10, -8, 3, -10, -12, -5, -13, -5, -9],
        [-20, 0, 0, -20, -10, -10, 0, -20, -10, -10, -13, -7, -7, -13, -10],
        [5, -20, -20, 5, -8, -8, -20, 5, -8, -8, -3, -12, -12, -3, -8],
        [5, 0, -20, -20, -10, 3, -10, -8, -20, -8, -5, -13, -5, -12, -9],
        [-20, 0, -20, 5, 3, -10, -10, -8, -8, -20, -5, -5, -13, -12, -9],
        [-3, -20, -7, -3, -12, -12, -13, -3, -5, -5, -9, -10, -10, -4, -8],
        [-3, -7, -7, -20, -13, -5, -7, -12, -13, -5, -10, -11, -6, -10, -9],
        [-20, -7, -7, -3, -5, -13, -7, -12, -5, -13, -10, -6, -11, -10, -9],
        [-3, -7, -20, -3, -5, -5, -13, -3, -12, -12, -4, -10, -10, -9, -8],
        [-8, -10, -10, -8, -9, -9, -10, -8, -9, -9, -8, -9, -9, -8, -9],
    ],
)
