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
    DEFAULT_MIN_OVERLAP,
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
    min_overlap: int = DEFAULT_MIN_OVERLAP

    def __post_init__(self) -> None:
        """Validate replisome configuration.

        Raises
            TypeError: If the primer is not a valid DNA primer.
        """
        if self.primer.dna_type != DNAType.PRIMER:
            raise TypeError("A target sequence had been used as a primer.")

        self.target = (
            self.target.pad(len(self.primer) - self.min_overlap).upper().reverse()
        )

        self.primer = self.primer.upper().reverse()
        self.__target_index_limit = slice(
            0, len(self.target) - len(self.primer) + self.min_overlap
        )
        self.__max_primability = 0
        self.__max_stability = 0
        self.__max_quality = 0

    def replicon_slice(self, k: int) -> slice:
        """The slice of the target DNA that is being simulated."""
        if k > self.target_index_limit.stop or k < 0:
            raise IndexError(
                f"Requested index {k} is out of range. (max: {self.target_index_limit.stop})"
            )
        return slice(k, k + len(self.primer))

    def calc_primability(self, k: int) -> float:
        """Calculate the primability of the primer."""
        return k

    @property
    def target_index_limit(self) -> slice:
        """Return the index limit of the target DNA."""
        return self.__target_index_limit

    @property
    def max_primability(self) -> float:
        """Calculate the primability of the primer."""
        return self.__max_primability

    @property
    def max_stability(self) -> float:
        """Calculate the stability of the primer."""
        return self.__max_stability

    @property
    def max_quality(self) -> float:
        """Calculate the quality of the primer."""
        return self.__max_quality
