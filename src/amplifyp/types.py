# -*- coding: utf-8 -*-
"""Amplify P - Types."""
from typing import Iterator


class NucleotidePairwiseWeightTbl:  # pylint: disable=too-few-public-methods
    """Nucleotide Pairwise Weight Table."""

    def __init__(self) -> None:
        """Construct a Nucleotide Pairwise Weight Table."""


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
