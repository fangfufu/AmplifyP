# -*- coding: utf-8 -*-
"""Amplify P - Replisome related."""

from functools import cache

from .dna import DNA, DNAType
from .origin import Origin
from .settings import (
    ReplicationConfig,
    DEFAULT_REPLICATION_CONFIG,
)


class Replisome:
    """A class representing a replisome.

    A replisome consists of a target DNA sequence and a primer DNA sequence. The
    replisome class represents the replisome complex - the DNA polymerase with
    the primer moves along the target sequence to find a good binding site.

    Attributes:
        target (DNA): The DNA sequence being replicated.
        primer (DNA): The DNA sequence that serves as a starting point for DNA
            synthesis.
        replication_config (ReplisomeConfig): The configuration of the replisome.

    Methods:
        origin_slice(k: int) -> slice: Returns a slice object that represents
            the range of the origin starting at index k.
        range_limit() -> slice: Return the index limit of the target DNA.
    """

    def __init__(
        self,
        target: DNA,
        primer: DNA,
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
            self.__target = target.circular_pad().upper().reverse()
        elif target.dna_type == DNAType.LINEAR:
            self.__target = target.upper().reverse()
        else:
            raise TypeError("Invalid DNA type for target sequence.")

        self.__primer = primer.upper().reverse()
        self.__replication_config = replication_config

        self.__range = range(0, len(self.target) - len(self.primer))

    @property
    def target(self) -> DNA:
        """Return the processed target DNA sequence."""
        return self.__target

    @property
    def primer(self) -> DNA:
        """Return the processed primer DNA sequence."""
        return self.__primer

    @property
    def replication_config(self) -> ReplicationConfig:
        """Return the configuration of the replisome."""
        return self.__replication_config

    @property
    def range(self) -> range:
        """Return the range limit of the target DNA."""
        return self.__range

    def __origin_slice(self, k: int) -> slice:
        """Returns a the origin slice object.

        Args:
            k (int): The starting index of the slice.

        Returns:
            slice: A slice object that represents the range of the target for
            a origin starting at index k.

        Raises:
            IndexError: If the requested index is out of range.
        """
        if k > self.range.stop or k < 0:
            raise IndexError(
                f"Requested index {k} is out of range. (max: {self.range.stop})"
            )
        return slice(k, k + len(self.primer))

    # WARNING: This might actually cause memory leak. We need to look into this
    # later.
    @cache  # pylint: disable=method-cache-max-size-none
    def origin(self, k: int) -> Origin:
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
            self.target[self.__origin_slice(k)].sequence,
            self.primer.sequence,
            self.replication_config,
        )
