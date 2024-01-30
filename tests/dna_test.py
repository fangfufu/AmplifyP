# -*- coding: utf-8 -*-
"""Tests related to the DNA class."""

from amplifyp.dna import DNA, DNAType, DNADirection, Primer


def test_dna_init() -> None:
    """Test the initialization of a DNA object with a given sequence."""
    dna = DNA("ATCG")
    assert dna.sequence == "ATCG"
    assert dna.type == DNAType.LINEAR
    assert dna.name == "ATCG"
    assert dna.direction == DNADirection.FWD


def test_dna_name_setter() -> None:
    """Test that the name of a DNA object can be set correctly."""
    dna = DNA("ATCG")
    dna.name = "test"
    assert dna.name == "test"


def test_dna_lower() -> None:
    """Test that the `lower` method of the `DNA` class."""
    dna = DNA("ATCG")
    assert dna.lower().sequence == "atcg"


def test_dna_upper() -> None:
    """Test the `upper` method of the `DNA` class."""
    dna = DNA("atcg")
    assert dna.upper().sequence == "ATCG"


def test_dna_complement() -> None:
    """Test the complement method of the DNA class."""
    dna = DNA("ATCG")
    assert dna.complement().sequence == "TAGC"


def test_dna_reverse() -> None:
    """Test the reverse method of the DNA class."""
    dna = DNA("ATCG")
    assert dna.reverse().sequence == "GCTA"


def test_dna_eq() -> None:
    """Test that the equality operator works for DNA objects."""
    dna1 = DNA("ATCG")
    dna2 = DNA("ATCG")
    assert dna1 == dna2


def test_dna_is_complement_of() -> None:
    """Test the is_complement_of method of the DNA class."""
    dna1 = DNA("ATCG")
    dna2 = DNA("TAGC", direction=DNADirection.REV)
    assert dna1.is_complement_of(dna2)


def test_dna_len() -> None:
    """Test the length of a DNA sequence."""
    dna = DNA("ATCG")
    assert len(dna) == 4


def test_dna_pad() -> None:
    """Test the padding method of the DNA class."""
    dna = DNA("ATCG", dna_type=DNAType.CIRCULAR)
    assert dna.pad(2).sequence == "CGATCG"
    dna = DNA("ATCG")
    assert dna.pad(2).sequence == "--ATCG"


def test_dna_getitem() -> None:
    """Test the __getitem__ method of the DNA class."""
    dna = DNA("ATCG")
    assert dna[1:3].sequence == "TC"


def test_dna_str() -> None:
    """Test the string representation of the DNA class."""
    dna = DNA("ATCG")
    assert str(dna) == "DNA: ATCG, LINEAR, FWD"


def test_primer_init() -> None:
    """Test the initialization of a Primer object."""
    primer = Primer("ATCG")
    assert primer.sequence == "ATCG"
    assert primer.type == DNAType.PRIMER
    assert primer.direction == DNADirection.FWD
    assert primer.name == "ATCG"
