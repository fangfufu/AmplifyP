# -*- coding: utf-8 -*-
"""Amplify P - DNA."""
from . import nucleotides


class DNA:
    """DNA sequence container."""

    def __init__(self, sequence: str, primer: bool = False) -> None:
        """Construct a DNA sequence."""
        check_str = nucleotides.PRIMER if primer else nucleotides.TARGET
        if not set(sequence.upper()) <= set(check_str):
            raise ValueError("DNA sequence contains invalid characters.")
        self._sequence = sequence

    def __str__(self) -> str:
        """Return the string representation of a DNA sequence."""
        return self.sequence

    @property
    def sequence(self) -> str:
        """Return the DNA sequence."""
        return self._sequence

    @property
    def reverse(self) -> str:
        """Return the reverse complement of the DNA sequence."""
        return self._sequence[::-1]

    @property
    def lower(self) -> str:
        """Return the DNA sequence in lower case."""
        return self._sequence.lower()

    @property
    def upper(self) -> str:
        """Return the DNA sequence in upper case."""
        return self._sequence.upper()

    @property
    def complement(self) -> str:
        """Return the complement of the DNA sequence.

        Note that the complement of non-ACGT bases are undefined.
        """
        return self._sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))
