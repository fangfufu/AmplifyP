# -*- coding: utf-8 -*-
"""Tests for replisome.py."""
from amplifyp.replisome import DNA


def test_dna_invalid_characters() -> None:
    """Test if DNA sequence contains invalid characters."""
    try:
        DNA("ACGTX")
    except ValueError as err:
        assert str(err) == "DNA sequence contains invalid characters."


def test_dna_str_representation() -> None:
    """Test the string representation of a DNA sequence."""
    dna = DNA("ACGT")
    assert str(dna) == "ACGT"


def test_dna_complement() -> None:
    """Test the complement of a DNA sequence."""
    dna = DNA("tACGTacgta")
    assert dna.complement == "aTGCAtgcat"


def test_dna_lower_case() -> None:
    """Test the lower case of a DNA sequence."""
    dna = DNA("ACGT")
    assert dna.lower == "acgt"


def test_dna_upper_case() -> None:
    """Test the upper case of a DNA sequence."""
    dna = DNA("acgt")
    assert dna.upper == "ACGT"


def test_dna_sequence_reverse() -> None:
    """Test the reverse of a DNA sequence."""
    dna = DNA("AACGTTA")
    assert dna.reverse == "ATTGCAA"
