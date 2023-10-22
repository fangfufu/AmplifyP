# -*- coding: utf-8 -*-
"""Amplify P - Replisome related."""


from .dna import DNA, DNAType, DNADirection
from .replication_origin import ReplicationOrigin
from .settings import (
    ReplicationConfig,
    DEFAULT_REPLICATION_CONFIG,
)


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
        self.__target_fwd: DNA = target.upper()
        if target.type == DNAType.CIRCULAR:
            self.__target_fwd = self.__target_fwd.circular_pad()

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

    def origin(self, k: int, primer: DNA, direction: DNADirection) -> ReplicationOrigin:
        """Returns an Origin object representing the origin of replication."""
        return ReplicationOrigin(
            # WARNING: The reveral here might cause performance issues
            self.target(direction).sequence[self.__origin_slice(k, primer)][::-1],
            primer.sequence[::-1],
            self.replication_config,
        )
