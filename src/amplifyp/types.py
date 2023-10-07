# -*- coding: utf-8 -*-
"""Amplify P - Types."""
from typing import Iterator


class NucleotidePairwiseWeightTbl:
    """Nucleotide Pairwise Weight Table."""

    def __init__(self, row: str, column: str, init_weight: float) -> None:
        """Construct a Nucleotide Pairwise Weight Table."""
        self._weight = dict.fromkeys([(x, y) for x in row for y in column], init_weight)

    def __getitem__(self, key: tuple[str, str]) -> float:
        """Return the weight of at certain nucleotide pair."""
        return self._weight[key]

    def __setitem__(self, key: tuple[str, str], value: float) -> None:
        """Set the weight at a certain nucleotide pair."""
        self._weight[key] = value

    def __len__(self) -> int:
        """Return the size of the Run-length Weight table."""
        return len(self._weight)

    def __str__(self) -> str:
        """Return the string representation of the table."""
        return str(self._weight)

    def __repr__(self) -> str:
        """Return the string representation of the table."""
        return str(self._weight)

    def __iter__(self) -> Iterator[tuple[str, str]]:
        """Return the iterator of the table."""
        return iter(self._weight)


class RunLengthWeightTbl:
    """Run-length Weight Table Class."""

    def __init__(self, size: int, init_weight: float) -> None:
        """Construct a Run-length Weight Table."""
        self._weight = [init_weight] * size

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
