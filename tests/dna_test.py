# -*- coding: utf-8 -*-
# pylint: disable=cyclic-import
"""Tests for types.py."""

import pytest

from amplifyp.dna import DNA, DNAType


def test_dna() -> None:
    """Test DNA class."""
    # Test the creation of a DNA with invalid symbols.
    with pytest.raises(ValueError):
        DNA("RYKMSW")

    # Test the creation of a DNA primer.
    seq = "RYKMSW"
    dna = DNA("RYKMSW", dna_type=DNAType.PRIMER)
    assert dna.sequence == seq

    # Test the complement of a DNA sequence.
    dna = DNA("ATGC")
    assert dna.complement().sequence == "GCAT"
    assert dna.complement() == dna

    # Test the equality of a DNA sequence.
    assert DNA("ATGC") == DNA("ATGC")
    assert DNA("ATGC") == DNA("GCAT")
    assert DNA("tACGTacgta") != DNA("AAAAAAAAA")

    # Test the lower case of a DNA sequence.
    dna = DNA("ACGT")
    assert dna.lower().sequence == "acgt"

    # Test the upper case of a DNA sequence.
    dna = DNA("acgt")
    assert dna.upper().sequence == "ACGT"

    # Test the concatenation of two DNA sequences.
    dna_1 = DNA("ACGT")
    dna_2 = DNA("TGCA")
    assert (dna_1 + dna_2).sequence == "ACGTTGCA"
