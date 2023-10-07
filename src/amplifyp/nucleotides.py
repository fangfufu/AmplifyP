# -*- coding: utf-8 -*-
"""Amplify P - Nucleotide constants."""
from typing import Final

# These are the IUPAC degenerate base symbols.

SINGLE: Final[str] = "GATC"
DOUBLE: Final[str] = "MRWSYK"
TRIPLE: Final[str] = "VHDB"
WILDCARD: Final[str] = "N"

PRIMER: Final[str] = SINGLE + DOUBLE + TRIPLE + WILDCARD
TARGET: Final[str] = SINGLE + WILDCARD
