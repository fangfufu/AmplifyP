# -*- coding: utf-8 -*-
# pylint: disable=cyclic-import
"""Tests for types.py."""

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
    """Test the pading method of the DNA class."""
    dna = DNA("ATCG", dna_type=DNAType.CIRCULAR)
    assert dna.pad(2).sequence == "CGATCG"
    dna = DNA("ATCG")
    assert dna.pad(2).sequence == "--ATCG"


def test_dna_getitem() -> None:
    """Test the __getitem__ method of the DNA class."""
    dna = DNA("ATCG")
    assert dna[1:3].sequence == "TC"


def test_dna_str() -> None:
    """Test the string representation of a DNA object."""
    dna = DNA("ATCG")
    assert str(dna) == "DNA: ATCG"


def test_dna_repr() -> None:
    """Test the __repr__ method of the DNA class."""
    dna = DNA("ATCG")
    assert repr(dna) == "DNA: Name: ATCG, DNAType.LINEAR, DNADir.FWD, Seq: ATCG"


def test_primer_init() -> None:
    """Test the initialization of a Primer object."""
    primer = Primer("ATCG")
    assert primer.sequence == "ATCG"
    assert primer.type == DNAType.PRIMER
    assert primer.direction == DNADirection.FWD
    assert primer.name == "ATCG"


def test_primer_index() -> None:
    """Test the index method of the Primer class."""
    dna = DNA("ATCG")
    primer = Primer("AT")
    primer.index.append(dna, DNADirection.FWD, 0)
    assert primer.index[dna, DNADirection.FWD] == [0]


def test_primer_index_clear() -> None:
    """Test the index_clear method."""
    dna = DNA("ATCG")
    primer = Primer("AT")
    primer.index.append(dna, DNADirection.FWD, 0)
    primer.index.clear(dna, DNADirection.FWD)
    assert not primer.index[dna, DNADirection.FWD]


def test_primer_index_remove() -> None:
    """Test the functionality of the Primer.index_remove() method."""
    dna = DNA("ATCG")
    primer = Primer("AT")
    primer.index.append(dna, DNADirection.FWD, 0)
    primer.index.remove(dna, DNADirection.FWD, 0)
    assert not primer.index[dna, DNADirection.FWD]


def test_primer_index_clear_all() -> None:
    """Test the functionality of the Primer.index.clear_all() method."""
    dna1 = DNA("ATCG")
    dna2 = DNA("ATCG")
    primer = Primer("AT")
    primer.index.append(dna1, DNADirection.FWD, 0)
    primer.index.append(dna2, DNADirection.FWD, 0)
    primer.index.clear_all()
    assert not primer.index[dna1, DNADirection.FWD]
    assert not primer.index[dna2, DNADirection.FWD]
