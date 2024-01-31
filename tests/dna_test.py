# -*- coding: utf-8 -*-
"""Tests related to the DNA class."""

import pytest

from amplifyp.dna import DNA, DNAType, DNADirection, Primer


def test_dna_init() -> None:
    """Test the initialization of a DNA object with a given sequence."""
    dna = DNA("ATCG")
    assert dna.seq == "ATCG"
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
    assert dna.lower().seq == "atcg"


def test_dna_upper() -> None:
    """Test the `upper` method of the `DNA` class."""
    dna = DNA("atcg")
    assert dna.upper().seq == "ATCG"


def test_dna_complement() -> None:
    """Test the complement method of the DNA class."""
    dna = DNA("ATCG")
    assert dna.complement().seq == "TAGC"


def test_dna_reverse() -> None:
    """Test the reverse method of the DNA class."""
    dna = DNA("ATCG")
    assert dna.reverse().seq == "GCTA"


def test_dna_eq() -> None:
    """Test that the equality operator works for DNA objects."""
    dna1 = DNA("ATCG")
    dna2 = DNA("ATCG")
    assert dna1 == dna2
    assert dna1 != ""


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
    assert dna.pad(2).seq == "CGATCG"
    dna = DNA("ATCG")
    assert dna.pad(2).seq == "--ATCG"
    with pytest.raises(TypeError):
        primer = Primer("ATCG")
        primer.pad(3)


def test_dna_getitem() -> None:
    """Test the __getitem__ method of the DNA class."""
    dna = DNA("ATCG")
    assert dna[1:3].seq == "TC"


def test_dna_str() -> None:
    """Test the string representation of the DNA class."""
    dna = DNA("ATCG")
    assert str(dna) == "DNA: ATCG, LINEAR, FWD"


def test_primer_init() -> None:
    """Test the initialization of a Primer object."""
    primer = Primer("ATCG")
    assert primer.seq == "ATCG"
    assert primer.type == DNAType.PRIMER
    assert primer.direction == DNADirection.FWD
    assert primer.name == "ATCG"


def test_dna_invalid_type() -> None:
    """Test DNA initialised with invalid type"""
    with pytest.raises(TypeError):
        DNA("A", 4)


def test_dna_invalid_char() -> None:
    """Test DNA initialised with invalid characters"""
    with pytest.raises(ValueError):
        DNA("L")


def test_dna_rotation() -> None:
    """Test DNA rotation."""
    a = DNA("AAAGG", DNAType.CIRCULAR)
    b = DNA("GGAAA", DNAType.CIRCULAR)
    c = DNA("G")
    assert a.rot(2) == b

    with pytest.raises(TypeError):
        c.rot(1)


def test_dna_hash() -> None:
    """Test that the DNA hash function generates a hash"""
    assert hash(DNA("A"))
