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
    dna = DNA("RYKMSW", dna_type=DNAType.PRIMER)
    assert dna.sequence == "RYKMSW"

    # Test the complement of a DNA sequence.
    dna = DNA("ACGTMKRYBDHVacgtmkrybdhv", dna_type=DNAType.PRIMER)
    assert dna.complement().sequence == "bdhvrymkacgtBDHVRYMKACGT"

    # Test the equality of a DNA sequence.
    assert DNA("ATGC") == DNA("ATGC")
    assert DNA("ATGC").is_complement_of(DNA("GCAT"))
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

    # Test the pad method.
    dna = DNA("ACGT", dna_type=DNAType.CIRCULAR)
    assert dna.circular_pad().sequence == "ACGTACGTACGT"
    dna = DNA("ACGT", dna_type=DNAType.LINEAR)
    assert dna.circular_pad() == NotImplemented

    # Test the length of a DNA sequence.
    dna = DNA("ACGT")
    assert len(dna) == 4

    # Test DNA sequencce setter
    dna = DNA("ACGT")
    assert dna[0].sequence == "A"
