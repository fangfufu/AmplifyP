# -*- coding: utf-8 -*-
"""Amplify P - settings related.

Don't ask me how the default values were chosen, I shamelessly copied them from
https://github.com/wrengels/Amplify4/blob/master/Amplify4/Amplify4/TargDelegate.swift
"""

from dataclasses import dataclass
from typing import Final, List, Tuple, Dict
from .dna import Nucleotides


class LengthWiseWeightTbl:
    """A class representing a run-length weight table."""

    def __init__(
        self,
        default_weight: float = 0,
        overrides: Dict[int, float] | None = None,
    ) -> None:
        """Initialize a Run-length Weight Table."""
        if overrides is None:
            overrides = {}
        self.__default_weight = default_weight
        self.__overrides = overrides

    def __getitem__(self, key: int) -> float:
        """Return the weight of at certain run-length."""
        if key in self.__overrides:
            return self.__overrides[key]
        return self.__default_weight

    def __setitem__(self, key: int, value: float) -> None:
        """Set the weight at a certain run-length."""
        self.__overrides[key] = value


class BasePairWeightsTbl:
    """Nucleotide Pairwise Weight Table."""

    def __init__(self, row: str, col: str, weight: List[List[float]]) -> None:
        """Construct a Nucleotide Pairwise Weight Table.

        Args:
            row (str): A string representing the row labels of the weight table.
            col (str): A string representing the column labels of the weight
                table.
            weight (List[List[float]]): A 2D list of floats representing the
                weights of each nucleotide pair.

        Raises:
            ValueError: If the length of the weight table does not match the
                length of the row or column labels.
        """
        self.__row = row
        self.__col = col
        self.__weight: Dict[Tuple[str, str], float] = {}
        self.__row_max: Dict[str, float] = {}

        # Expected row and column lengths
        exp_row_len = len(row) if Nucleotides.GAP not in row else len(row) - 1
        exp_col_len = len(col) if Nucleotides.GAP not in col else len(col) - 1

        if len(weight) != exp_row_len:
            raise ValueError(
                "BasePairWeightsTbl: row length mismatch at initialisation."
            )

        for i, row_val in enumerate(self.__row):
            if row_val != Nucleotides.GAP:
                # We never put the gap symbol in the table, hence the -1.
                if len(weight[i]) != exp_col_len:
                    raise ValueError(
                        "BasePairWeightsTbl: column length mismatch at initialisation."
                    )
                self.__row_max[row_val] = max(weight[i])
            for j, col_val in enumerate(self.__col):
                if Nucleotides.GAP in [row_val, col_val]:
                    self.__weight[row_val, col_val] = 0
                else:
                    self.__weight[row_val, col_val] = weight[i][j]

    def row(self) -> str:
        """Return the row nucleotides."""
        return self.__row[:-1]

    def column(self) -> str:
        """Return the column nucleotides."""
        return self.__col[:-1]

    def row_max(self, row: str) -> float:
        """Return the maximum weight of a row."""
        return self.__row_max[row]

    def __getitem__(self, key: tuple[str, str]) -> float:
        """Return the weight of at certain nucleotide pair."""
        i, j = key
        i = i.upper()
        j = j.upper()
        return self.__weight[i, j]

    def __setitem__(self, key: tuple[str, str], value: float) -> None:
        """Set the weight at a certain nucleotide pair."""
        i, j = key
        i = i.upper()
        j = j.upper()
        self.__weight[key] = value

    def __len__(self) -> int:
        """Return the size of the Run-length Weight table."""
        return len(self.row()) * len(self.column())

    def __str__(self) -> str:
        """Return the string representation of the table."""
        return str(self.__weight)


DEFAULT_MATCH_WEIGHTS: Final[LengthWiseWeightTbl] = LengthWiseWeightTbl(
    default_weight=1,
    overrides={
        0: 30,
        1: 20,
        2: 10,
        3: 10,
        4: 9,
        5: 9,
        6: 8,
        7: 7,
        8: 6,
        9: 5,
        10: 5,
        11: 4,
        12: 3,
        13: 2,
        14: 1,
    },
)
DEFAULT_RUN_WEIGHTS: Final[LengthWiseWeightTbl] = LengthWiseWeightTbl(
    default_weight=186,
    overrides={0: 100, 1: 150, 2: 175, 3: 182, 4: 186},
)

DEFAULT_BASE_PAIR_WEIGHTS: Final[BasePairWeightsTbl] = BasePairWeightsTbl(
    row=Nucleotides.PRIMER,
    col=Nucleotides.LINEAR,
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
    col=Nucleotides.PRIMER,
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

# Default threshold for primability.
DEFAULT_PRIMABILITY_CUTOFF: Final[float] = 0.8

# Default threshold for stability.
DEFAULT_STABILITY_CUTOFF: Final[float] = 0.4


@dataclass(slots=True)
class Settings:
    """Configuration class for replication settings.

    Attributes:
        base_pair_scores (BasePairWeightsTbl): Table of weights for each base pair.
        match_weight (LengthWiseWeightTbl): Table of weights for each match.
        run_weight (LengthWiseWeightTbl): Table of weights for each run.
        primability_cutoff (float): Cutoff value for primability.
        stability_cutoff (float): Cutoff value for stability.
    """

    base_pair_scores: BasePairWeightsTbl = DEFAULT_BASE_PAIR_WEIGHTS
    match_weight: LengthWiseWeightTbl = DEFAULT_MATCH_WEIGHTS
    run_weights: LengthWiseWeightTbl = DEFAULT_RUN_WEIGHTS
    primability_cutoff: float = DEFAULT_PRIMABILITY_CUTOFF
    stability_cutoff: float = DEFAULT_STABILITY_CUTOFF


DEFAULT_SETTINGS: Final[Settings] = Settings()
