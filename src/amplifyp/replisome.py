# -*- coding: utf-8 -*-
"""Amplify P - Data types."""
from dataclasses import dataclass

from .dna import DNA, DNAType
from .settings import (
    ReplisomeConfig,
    DEFAULT_REPLISOME_CONFIG,
)


@dataclass
class Replicon:
    """A class representing the region of DNA being replicated.

    Attributes:
        target (str): The target DNA sequence.
        primer (str): The primer sequence.
        replisome_config (ReplisomeConfig): The configuration for the replisome.
            Defaults to DEFAULT_REPLISOME_CONFIG.

    Raises:
        ValueError: If the length of the target is not equal to the length of the primer.
    """

    target: str
    primer: str
    replisome_config: ReplisomeConfig = DEFAULT_REPLISOME_CONFIG

    def __post_init__(self) -> None:
        """Validates that the length of the target and primer are equal."""
        if len(self.target) != len(self.primer):
            raise ValueError("The target has to have the same length as the primer.")


@dataclass
class Replisome:
    """A class representing a replisome.

    A replisome consists of a target DNA sequence and a primer DNA sequence. The
    replisome class represents the replisome complex - the DNA polymerase with
    the primer moves along the target sequence to find a good binding site.

    Attributes:
        target (DNA): The DNA sequence being replicated.
        primer (DNA): The DNA sequence that serves as a starting point for DNA
            synthesis.
        replisome_config (ReplisomeConfig): The configuration of the replisome.

    Methods:
        replicon_slice(k: int) -> slice: Returns a slice object that represents
            the range of the replicon starting at index k.
        range_limit() -> slice: Return the index limit of the target DNA.
    """

    target: DNA
    primer: DNA
    replisome_config: ReplisomeConfig = DEFAULT_REPLISOME_CONFIG

    def __post_init__(self) -> None:
        """Validate replisome configuration.

        Raises
            TypeError: If the primer is not a valid DNA primer.
        """
        if self.primer.dna_type != DNAType.PRIMER:
            raise TypeError("A target sequence had been used as a primer.")

        self.target = (
            self.target.pad(len(self.primer) - self.replisome_config.min_overlap)
            .upper()
            .reverse()
        )

        self.primer = self.primer.upper().reverse()
        self.__range_limit = range(
            0, len(self.target) - len(self.primer) + self.replisome_config.min_overlap
        )
        self.__max_primability = 0
        self.__max_stability = 0
        self.__max_quality = 0

    @property
    def range_limit(self) -> range:
        """Return the index limit of the target DNA."""
        return self.__range_limit

    def replicon_slice(self, k: int) -> slice:
        """Returns a the replicon slice object.

        Args:
            k (int): The starting index of the slice.

        Returns:
            slice: A slice object that represents the range of the target for
            a replicon starting at index k.

        Raises:
            IndexError: If the requested index is out of range.
        """
        if k > self.range_limit.stop or k < 0:
            raise IndexError(
                f"Requested index {k} is out of range. (max: {self.range_limit.stop})"
            )
        return slice(k, k + len(self.primer))

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
