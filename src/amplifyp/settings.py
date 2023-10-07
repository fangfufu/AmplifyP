# -*- coding: utf-8 -*-
"""Amplify P - Settings."""

from typing import Final

from . import nucleotides
from .types import LengthWiseWeightTbl, NucleotidePairwiseWeightTbl

DEFAULT_MATCH_WEIGHTS_LEN: Final[int] = 200
DEFAULT_RUN_WEIGHTS_LEN: Final[int] = 200

base_pair_scores: NucleotidePairwiseWeightTbl = NucleotidePairwiseWeightTbl(
    row=nucleotides.PRIMER, column=nucleotides.TARGET
)

match_weights: LengthWiseWeightTbl = LengthWiseWeightTbl(DEFAULT_MATCH_WEIGHTS_LEN)

run_weights: LengthWiseWeightTbl = LengthWiseWeightTbl(DEFAULT_RUN_WEIGHTS_LEN)

primer_dimer_weights: NucleotidePairwiseWeightTbl = NucleotidePairwiseWeightTbl(
    row=nucleotides.PRIMER, column=nucleotides.PRIMER
)
