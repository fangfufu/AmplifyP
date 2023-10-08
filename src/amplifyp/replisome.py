# -*- coding: utf-8 -*-
"""Amplify P - Data types."""
from dataclasses import dataclass

from .dna import DNA, DNAType
from .settings import (
    BasePairWeights,
    LengthWiseWeightTbl,
    DEFAULT_BASE_PAIR_WEIGHTS,
    DEFAULT_MATCH_WEIGHTS,
    DEFAULT_RUN_WEIGHTS,
)


@dataclass(frozen=True, slots=True)
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
    base_pair_scores: BasePairWeights = DEFAULT_BASE_PAIR_WEIGHTS
    match_weight: LengthWiseWeightTbl = DEFAULT_MATCH_WEIGHTS
    run_weight: LengthWiseWeightTbl = DEFAULT_RUN_WEIGHTS

    def __post_init__(self) -> None:
        """Validate replisome configuration.

        Raises
            TypeError: If the primer is not a valid DNA primer.
        """
        if self.primer.dna_type != DNAType.PRIMER:
            raise TypeError("A target sequence had been used as a primer.")

    def primability(self, k: int) -> float:
        """Calculate the primability of the primer at target position k."""
        return k

    def stability(self, k: int) -> float:
        """Calculate the stability of the primer at target position k."""
        return k

    def quality(self, k: int) -> float:
        """Calculate the quality of the primer at target position k."""
        return k
