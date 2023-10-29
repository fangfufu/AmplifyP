# -*- coding: utf-8 -*-
"""Amplify P - Replication related."""

from dataclasses import dataclass
from functools import cached_property
from typing import List

from .dna import DNA, DNADirection
from .settings import (
    ReplicationConfig,
    LengthWiseWeightTbl,
    BasePairWeightsTbl,
    DEFAULT_REPLICATION_CONFIG,
)


@dataclass(frozen=True)
class ReplicationOrigin:
    """A class representing the origin of replication.

    Attributes:
        target (str): The target DNA sequence as a string, in 3'-5' orientation.
        primer (str): The primer sequence as a string, in 3'-5' orientation.
        replication_config (ReplisomeConfig): The configuration for the replisome.
            Defaults to DEFAULT_REPLISOME_CONFIG.

    Raises:
        ValueError: If the length of the target is not equal to the length of the primer.
    """

    target: str
    primer: str
    replication_config: ReplicationConfig

    def __post_init__(self) -> None:
        """Validates that the length of the target and primer are equal."""
        if len(self.target) != len(self.primer):
            raise ValueError("The target has to have the same length as the primer.")

    def primability(self) -> float:
        """Returns the primability of the origin.

        Returns:
            float: The primability of the origin.
        """
        m: LengthWiseWeightTbl = self.replication_config.match_weight
        S: BasePairWeightsTbl = (  # pylint: disable=invalid-name
            self.replication_config.base_pair_scores
        )
        numerator: float = 0
        denominator: float = 0
        for k, (i, j) in enumerate(zip(self.primer, self.target)):
            numerator += m[k] * S[i, j]
            denominator += m[k] * S.row_max(i)
        score = numerator / denominator
        return score

    def stability(self) -> float:
        """Returns the stability of the origin.

        Returns:
            float: The stability of the origin.
        """
        r = self.replication_config.run_weight
        S = self.replication_config.base_pair_scores  # pylint: disable=invalid-name
        numerator: float = 0
        denominator: float = 0
        Rn: float = 0  # pylint: disable=invalid-name
        for k, (i, j) in enumerate(zip(self.primer, self.target)):
            numerator += r[k] * S[i, j]
            denominator += S.row_max(i)
            if r[k] > Rn:
                Rn = r[k]  # pylint: disable=invalid-name
            if S[i, j] <= 0:
                break
        score = numerator / (Rn * denominator)
        return score

    def quality(self) -> float:
        """Returns the quality of the origin.

        Returns:
            float: The quality of the origin.
        """
        cutoffs = (
            self.replication_config.primability_cutoff
            + self.replication_config.stability_cutoff
        )
        return (self.primability() + self.stability() - cutoffs) / (2 - cutoffs)


class ReplicationTarget:
    """
    A class representing a replisome.

    A replisome is a complex of proteins that carries out DNA replication. This
    class provides a representation of a replisome, including the target DNA
    sequence to be replicated, the configuration for DNA replication, and methods
    for accessing and manipulating the target sequence.
    """

    def __init__(
        self,
        target: DNA,
        replication_config: ReplicationConfig = DEFAULT_REPLICATION_CONFIG,
    ):
        """
        Initializes a new instance of the Replicator class.

        Args:
            target (DNA): The DNA sequence to be replicated.
            replication_config (ReplicationConfig, optional):
                The configuration for the replication process. Defaults to
                DEFAULT_REPLICATION_CONFIG.
        """
        self.target = target

        self.__replication_config: ReplicationConfig = replication_config

    @property
    def target(self) -> DNA:
        """Return the target DNA sequence."""
        return self.__target_fwd

    @target.setter
    def target(self, target: DNA) -> None:
        self.__target_fwd: DNA = target.upper()

    @cached_property
    def target_rev(self) -> DNA:
        """Return the reverse complement of the target DNA sequence."""
        return self.__target_fwd.complement()

    @property
    def replication_config(self) -> ReplicationConfig:
        """Return the configuration of the replisome."""
        return self.__replication_config

    @replication_config.setter
    def replication_config(self, config: ReplicationConfig) -> None:
        """Set the configuration of the replisome."""
        self.__replication_config = config

    def target_range(self, primer: DNA) -> range:
        """Return the range limit for the target DNA, given a primer sequence."""
        return range(0, len(self.target) - len(primer))

    def __origin_slice(self, k: int, primer: DNA) -> slice:
        """Returns a slice object for extracting the replication origin."""
        target_range = self.target_range(primer)
        if k > target_range.stop or k < 0:
            raise IndexError(
                f"Requested index {k} is out of range. (max: {target_range.stop})"
            )
        return slice(k, k + len(primer))

    def origin(self, k: int, primer: DNA, direction: DNADirection) -> ReplicationOrigin:
        """Returns an Origin object representing the origin of replication."""
        target: DNA = (
            self.target if direction == DNADirection.FORWARD else self.target_rev
        )
        return ReplicationOrigin(
            # WARNING: The reversal here might cause performance issues
            target.sequence[self.__origin_slice(k, primer)][::-1],
            primer.sequence[::-1],
            self.replication_config,
        )


@dataclass(frozen=True, slots=True)
class PrimerScore:
    """
    A class representing a primer's score.

    Attributes:
        primer (DNA): The DNA sequence of the primer.
        index (int): The index of the primer.
        primability (float): The primability score of the primer.
        stability (float): The stability score of the primer.
        quality (float): The quality score of the primer.
    """

    primer: DNA
    index: int
    primability: float
    stability: float
    quality: float


class Replisome:
    """A class representing a DNA replication complex, or replisome."""

    def __init__(
        self,
        target: DNA,
        primers: List[DNA],
        replication_config: ReplicationConfig = DEFAULT_REPLICATION_CONFIG,
    ):
        """
        Initializes a Replication object.

        Args:
            target (DNA): The DNA sequence to be replicated.
            primers (List[DNA]): A list of DNA sequences to be used as primers
                for replication.
            replication_config (ReplicationConfig, optional): The configuration
                for the replication process. Defaults to
                DEFAULT_REPLICATION_CONFIG.
        """
        self.__replication_target = ReplicationTarget(target, replication_config)
        self.__primers = primers

    @property
    def target(self) -> DNA:
        """Return the target DNA sequence."""
        return self.__replication_target.target

    @target.setter
    def target(self, value: DNA) -> None:
        """Set the target DNA sequence."""
        self.__replication_target.target = value

    @property
    def primers(self) -> List[DNA]:
        """Return the list of primers."""
        return self.__primers

    @primers.setter
    def primers(self, value: List[DNA]) -> None:
        """Set the list of primers."""
        self.__primers = value

    @property
    def replication_config(self) -> ReplicationConfig:
        """Return the configuration of the replisome."""
        return self.__replication_target.replication_config

    @replication_config.setter
    def replication_config(self, value: ReplicationConfig) -> None:
        """Set the configuration of the replisome."""
        self.__replication_target.replication_config = value
