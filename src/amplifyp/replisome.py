# -*- coding: utf-8 -*-
"""Amplify P - Data types."""
from dataclasses import dataclass

from .dna import DNA, DNAType
from .settings import (
    BasePairWeightsTbl,
    LengthWiseWeightTbl,
    DEFAULT_BASE_PAIR_WEIGHTS,
    DEFAULT_MATCH_WEIGHTS,
    DEFAULT_RUN_WEIGHTS,
)


@dataclass
class Replisome:  # pylint: disable=too-many-instance-attributes
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
    base_pair_scores: BasePairWeightsTbl = DEFAULT_BASE_PAIR_WEIGHTS
    match_weight: LengthWiseWeightTbl = DEFAULT_MATCH_WEIGHTS
    run_weight: LengthWiseWeightTbl = DEFAULT_RUN_WEIGHTS

    def __post_init__(self) -> None:
        """Validate replisome configuration.

        Raises
            TypeError: If the primer is not a valid DNA primer.
        """
        if self.primer.dna_type != DNAType.PRIMER:
            raise TypeError("A target sequence had been used as a primer.")

        self.target = self.target.pad(len(self.primer)).upper()
        self.primer = self.primer.upper()
        self.target_limits = range(0, len(self.target) - len(self.primer) - 1)
        self.__primability = 0
        self.__stability = 0
        self.__quality = 0

    def _calc_primability(self, k: int) -> float:
        """Calculate the primability of the primer."""
        return k

    @property
    def primability(self) -> float:
        """Calculate the primability of the primer."""
        return self.__primability

    @property
    def stability(self) -> float:
        """Calculate the stability of the primer."""
        return self.__stability

    @property
    def quality(self) -> float:
        """Calculate the quality of the primer."""
        return self.__quality
