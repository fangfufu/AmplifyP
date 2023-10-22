# -*- coding: utf-8 -*-
"""Amplify P - Replisome related"""


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
        if target.dna_type == DNAType.CIRCULAR:
            self.__target_fwd = target.circular_pad().upper()
        elif target.dna_type == DNAType.LINEAR:
            self.__target_fwd = target.upper().reverse()
        else:
            raise TypeError("Invalid DNA type for target sequence.")
        self.__target_rev: DNA = self.__target_fwd.complement().reverse()

        self.__replication_config: ReplicationConfig = replication_config

    def target(self, direction: DNADirection = DNADirection.FORWARD) -> DNA:
        """Return the processed target DNA sequence.

        Args:
            direction (DNADirection): The direction of the DNA sequence.
                Defaults to DNADirection.FORWARD.

        Returns:
            DNA: The processed target DNA sequence.
        """

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
        """Return the range limit of the target DNA."""
        return range(0, len(self.target()) - len(primer))

    def __origin_slice(self, k: int, primer: DNA) -> slice:
        """Returns a the origin slice object.

        Args:
            k (int): The starting index of the slice.

        Returns:
            slice: A slice object that represents the range of the target for
            a origin starting at index k.

        Raises:
            IndexError: If the requested index is out of range.
        """
        target_range = self.target_range(primer)
        if k > target_range.stop or k < 0:
            raise IndexError(
                f"Requested index {k} is out of range. (max: {target_range.stop})"
            )
        return slice(k, k + len(primer))

    def origin(self, k: int, primer: DNA) -> Origin:
        """Returns a origin object.

        Args:
            k (int): The starting index of the origin.

        Returns:
            Origin: A origin object that represents the target and primer
            sequences for a origin starting at index k.

        Raises:
            IndexError: If the requested index is out of range.
        """
        return Origin(
            self.target()[self.__origin_slice(k, primer)].sequence,
            # WARNING: this *might* cause performance issues.
            primer.sequence[::-1],
            self.replication_config,
        )
