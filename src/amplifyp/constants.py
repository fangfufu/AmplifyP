# -*- coding: utf-8 -*-
"""Amplify P - Constants."""

from dataclasses import dataclass


@dataclass(frozen=True)
class NucleotideSymbolContainer:
    """Nucleotide symbol container."""

    single: str
    double: str
    triple: str
    wildcard: str


Nucleotide: NucleotideSymbolContainer = NucleotideSymbolContainer(
    "GATC", "MRWSYK", "VHDB", "N"
)
