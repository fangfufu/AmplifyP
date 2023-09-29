# -*- coding: utf-8 -*-
from typing import Final, Set


class DNA:
    """Contains a DNA sequence."""

    def __init__(self, sequence: str):
        """Constructor for DNA class.

        We sanitize the sequence by making sure it does not contain any
        invalid characters.
        """
        sequence = sequence.upper()
        valid_symbol: Final[Set[str]] = {
            "A",  # Single bases
            "C",
            "G",
            "T",
            "U",
            "R",  # Double bases
            "Y",
            "K",
            "M",
            "S",
            "W",
            "B",  # Triple codes
            "D",
            "H",
            "V",
            "N",  # Wildcard
        }
        if sequence not in valid_symbol:
            raise ValueError("Invalid DNA sequence.")

        self.sequence = sequence
