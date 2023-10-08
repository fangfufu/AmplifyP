# -*- coding: utf-8 -*-
"""Amplify P - Settings.

Don't ask me how these values were chosen, I copied them from
https://github.com/wrengels/Amplify4/blob/master/Amplify4/Amplify4/TargDelegate.swift
"""

from typing import Final
from .types import Nucleotides, LengthWiseWeightTbl, BasePairWeights


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

DEFAULT_BASE_PAIR_WEIGHTS: Final[BasePairWeights] = BasePairWeights(
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
DEFAULT_PRIMER_DIMER_WEIGHTS: Final[BasePairWeights] = BasePairWeights(
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
        [50, -20, -20, -10, 3, -10, -8, -20, -8, -5, -13, -5, -12, -9],
        [-20, 0, -20, 5, 3, -10, -10, -8, -8, -20, -5, -5, -13, -12, -9],
        [-3, -20, -7, -3, -12, -12, -13, -3, -5, -5, -9, -10, -10, -4, -8],
        [-3, -7, -7, -20, -13, -5, -7, -12, -13, -5, -10, -11, -6, -10, -9],
        [-20, -7, -7, -3, -5, -13, -7, -12, -5, -13, -10, -6, -11, -10, -9],
        [-3, -7, -20, -3, -5, -5, -13, -3, -12, -12, -4, -10, -10, -9, -8],
        [-8, -10, -10, -8, -9, -9, -10, -8, -9, -9, -8, -9, -9, -8, -9],
    ],
)
