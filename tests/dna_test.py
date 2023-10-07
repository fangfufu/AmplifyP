# -*- coding: utf-8 -*-
"""Tests for dna.py."""
from amplifyp.dna import DNA, Nucleotides


def test_nucleotide_is_primer_sequence() -> None:
    """Test function for is_primer_sequence()."""
    assert Nucleotides.is_primer_sequence("GATC")
    assert Nucleotides.is_primer_sequence("MRWSYK")
    assert Nucleotides.is_primer_sequence("VHDB")
    assert Nucleotides.is_primer_sequence("N")
    assert Nucleotides.is_primer_sequence("GATCMRWSYKVHDBN")
    assert not Nucleotides.is_primer_sequence("GATCX")


def test_nucleotide_is_target_sequence() -> None:
    """Test function for is_target_sequence()."""
    assert Nucleotides.is_target_sequence("GATC")
    assert not Nucleotides.is_target_sequence("MRWSYK")
    assert not Nucleotides.is_target_sequence("VHDB")
    assert Nucleotides.is_target_sequence("N")
    assert Nucleotides.is_target_sequence("GATCN")
    assert not Nucleotides.is_target_sequence("GATCX")


def test_dna_invalid_characters() -> None:
    """Test the creation of a DNA with invalid symbols."""
    try:
        DNA("RYKMSW")
    except ValueError as err:
        assert str(err) == "DNA sequence contains invalid characters."


def test_dna_primer_creation() -> None:
    """Test the creation of a DNA primer."""
    dna = DNA("RYKMSW", primer=True)
    assert dna.sequence == "RYKMSW"


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
