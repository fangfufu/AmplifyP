# -*- coding: utf-8 -*-
"""Amplify P - Settings.

Don't ask me how the values were chosen, I copied them from
https://github.com/wrengels/Amplify4/commit/6f9a88e30c352fc73990cc3360e2de34b4279a1d
"""

from typing import Final

from . import nucleotides
from .types import LengthWiseWeightTbl, NucleotidePairwiseWeightTbl

DEFAULT_MATCH_WEIGHTS_LEN: Final[int] = 200
DEFAULT_RUN_WEIGHTS_LEN: Final[int] = 200

INIT_MATCH_WEIGHTS: Final[float] = 1
INIT_RUN_WEIGHT: Final[float] = 186

match_weights: None | LengthWiseWeightTbl = None
run_weights: None | LengthWiseWeightTbl = None
base_pair_scores: NucleotidePairwiseWeightTbl | None = None
primer_dimer_weights: NucleotidePairwiseWeightTbl | None = None


def reset_match_weights() -> None:
    """Reset Match Weight to a default set of values."""
    global match_weights  # pylint: disable=global-statement
    match_weights = LengthWiseWeightTbl(DEFAULT_MATCH_WEIGHTS_LEN, 1)
    match_weights[0] = 30
    match_weights[1] = 20
    match_weights[2] = 10
    match_weights[3] = 10
    match_weights[4] = 9
    match_weights[5] = 9
    match_weights[6] = 8
    match_weights[7] = 7
    match_weights[8] = 6
    match_weights[9] = 5
    match_weights[10] = 5
    match_weights[11] = 4
    match_weights[12] = 3
    match_weights[13] = 2
    match_weights[14] = INIT_MATCH_WEIGHTS


def reset_run_weights() -> None:
    """Reset Run Weights to a default set of values."""
    global run_weights  # pylint: disable=global-statement
    run_weights = LengthWiseWeightTbl(DEFAULT_RUN_WEIGHTS_LEN)
    run_weights[0] = 100
    run_weights[1] = 150
    run_weights[2] = 175
    run_weights[3] = 182
    run_weights[3] = INIT_RUN_WEIGHT


def reset_base_pair_scores() -> None:
    """Reset Base Pair Scores."""
    # pylint: disable=too-many-statements
    global base_pair_scores  # pylint: disable=global-statement
    base_pair_scores = NucleotidePairwiseWeightTbl(
        row=nucleotides.PRIMER, column=nucleotides.TARGET
    )

    row: str = "G"
    base_pair_scores[(row, "G")] = 100
    base_pair_scores[(row, "A")] = 0
    base_pair_scores[(row, "T")] = 0
    base_pair_scores[(row, "C")] = 0
    base_pair_scores[(row, "N")] = 30

    row = "A"
    base_pair_scores[(row, "G")] = 0
    base_pair_scores[(row, "A")] = 100
    base_pair_scores[(row, "T")] = 0
    base_pair_scores[(row, "C")] = 0
    base_pair_scores[(row, "N")] = 30

    row = "T"
    base_pair_scores[(row, "G")] = 0
    base_pair_scores[(row, "A")] = 0
    base_pair_scores[(row, "T")] = 100
    base_pair_scores[(row, "C")] = 0
    base_pair_scores[(row, "N")] = 30

    row = "C"
    base_pair_scores[(row, "G")] = 0
    base_pair_scores[(row, "A")] = 0
    base_pair_scores[(row, "T")] = 0
    base_pair_scores[(row, "C")] = 100
    base_pair_scores[(row, "N")] = 30

    row = "M"
    base_pair_scores[(row, "G")] = 0
    base_pair_scores[(row, "A")] = 70
    base_pair_scores[(row, "T")] = 0
    base_pair_scores[(row, "C")] = 70
    base_pair_scores[(row, "N")] = 30

    row = "R"
    base_pair_scores[(row, "G")] = 70
    base_pair_scores[(row, "A")] = 70
    base_pair_scores[(row, "T")] = 0
    base_pair_scores[(row, "C")] = 0
    base_pair_scores[(row, "N")] = 30

    row = "W"
    base_pair_scores[(row, "G")] = 0
    base_pair_scores[(row, "A")] = 70
    base_pair_scores[(row, "T")] = 70
    base_pair_scores[(row, "C")] = 0
    base_pair_scores[(row, "N")] = 30

    row = "S"
    base_pair_scores[(row, "G")] = 70
    base_pair_scores[(row, "A")] = 0
    base_pair_scores[(row, "T")] = 0
    base_pair_scores[(row, "C")] = 70
    base_pair_scores[(row, "N")] = 30

    row = "Y"
    base_pair_scores[(row, "G")] = 0
    base_pair_scores[(row, "A")] = 0
    base_pair_scores[(row, "T")] = 70
    base_pair_scores[(row, "C")] = 70
    base_pair_scores[(row, "N")] = 30

    row = "K"
    base_pair_scores[(row, "G")] = 70
    base_pair_scores[(row, "A")] = 0
    base_pair_scores[(row, "T")] = 70
    base_pair_scores[(row, "C")] = 0
    base_pair_scores[(row, "N")] = 30

    row = "V"
    base_pair_scores[(row, "G")] = 50
    base_pair_scores[(row, "A")] = 50
    base_pair_scores[(row, "T")] = 0
    base_pair_scores[(row, "C")] = 50
    base_pair_scores[(row, "N")] = 30

    row = "H"
    base_pair_scores[(row, "G")] = 0
    base_pair_scores[(row, "A")] = 50
    base_pair_scores[(row, "T")] = 50
    base_pair_scores[(row, "C")] = 50
    base_pair_scores[(row, "N")] = 30

    row = "D"
    base_pair_scores[(row, "G")] = 50
    base_pair_scores[(row, "A")] = 50
    base_pair_scores[(row, "T")] = 50
    base_pair_scores[(row, "C")] = 0
    base_pair_scores[(row, "N")] = 30

    row = "B"
    base_pair_scores[(row, "G")] = 50
    base_pair_scores[(row, "A")] = 0
    base_pair_scores[(row, "T")] = 50
    base_pair_scores[(row, "C")] = 50
    base_pair_scores[(row, "N")] = 30

    row = "N"
    base_pair_scores[(row, "G")] = 30
    base_pair_scores[(row, "A")] = 30
    base_pair_scores[(row, "T")] = 30
    base_pair_scores[(row, "C")] = 30
    base_pair_scores[(row, "N")] = 30


def reset_primer_dimer_weights() -> None:
    """Reset Base Pair Scores for primer dimer evaluation."""
    global primer_dimer_weights  # pylint: disable=global-statement
    primer_dimer_weights = NucleotidePairwiseWeightTbl(
        row=nucleotides.PRIMER, column=nucleotides.PRIMER
    )
