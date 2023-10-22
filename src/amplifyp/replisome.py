# -*- coding: utf-8 -*-
"""Amplify P - Replisome related."""


from .dna import DNA, DNAType, DNADirection
from .origin import Origin
from .settings import (
    ReplicationConfig,
    DEFAULT_REPLICATION_CONFIG,
)


class Replisome:
    """A class representing a replisome."""

    def __init__(
        self,
        target: DNA,
        replication_config: ReplicationConfig = DEFAULT_REPLICATION_CONFIG,
    ):
        """Initializes a Replisome object.

        Args:
            target (DNA): The target DNA sequence to replicate.
            primer (DNA): The primer DNA sequence to use for replication.
            replication_config (ReplisomeConfig, optional): The configuration
                for the replisome. Defaults to DEFAULT_REPLISOME_CONFIG.

        Raises:
            TypeError: If the primer sequence is not a primer DNA sequence.
        """
        if target.type == DNAType.CIRCULAR:
            self.__target_fwd = target.circular_pad().upper()
        elif target.type == DNAType.LINEAR:
            self.__target_fwd = target.upper()
        else:
            raise TypeError("Invalid DNA type for target sequence.")
        self.__target_rev: DNA = self.__target_fwd.complement()

        self.__replication_config: ReplicationConfig = replication_config

    def target(self, direction: DNADirection) -> DNA:
        """Return the target DNA sequence."""
        return (
            self.__target_fwd
            if direction == DNADirection.FORWARD
            else self.__target_rev
        )

    @property
    def replication_config(self) -> ReplicationConfig:
        """Return the configuration of the replisome."""
        return self.__replication_config

    def target_range(self, primer: DNA) -> range:
        """Return the range limit for the target DNA, given a primer sequence."""
        return range(0, len(self.target(DNADirection.FORWARD)) - len(primer))

    def __origin_slice(self, k: int, primer: DNA) -> slice:
        """Returns a slice object for extracting the replication origin."""
        target_range = self.target_range(primer)
        if k > target_range.stop or k < 0:
            raise IndexError(
                f"Requested index {k} is out of range. (max: {target_range.stop})"
            )
        return slice(k, k + len(primer))

    def origin(self, k: int, primer: DNA, direction: DNADirection) -> Origin:
        """Returns an Origin object representing the origin of replication."""
        return Origin(
            # WARNING: The reveral here might cause performance issues
            self.target(direction).sequence[self.__origin_slice(k, primer)][::-1],
            primer.sequence[::-1],
            self.replication_config,
        )
