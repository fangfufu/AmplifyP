# -*- coding: utf-8 -*-
# pylint: disable=cyclic-import
"""Tests for types.py."""

from amplifyp.dna import DNA, DNAType, DNADirection


def test_dna_init() -> None:
    """
    Test the initialization of a DNA object with a given sequence.
    """
    dna = DNA("ATCG")
    assert dna.sequence == "ATCG"
    assert dna.type == DNAType.LINEAR
    assert dna.name == "ATCG"
    assert dna.direction == DNADirection.FORWARD


def test_dna_sequence_setter() -> None:
    """
    Test the sequence setter method of the DNA class.

    Creates a DNA object with the sequence "ATCG", sets the sequence to "CGTA",
    and checks that the sequence was correctly updated.
    """
    dna = DNA("ATCG")
    dna.sequence = "CGTA"
    assert dna.sequence == "CGTA"


def test_dna_type_setter() -> None:
    """
    Test that the type of a DNA sequence can be set using the type setter method.
    """
    dna = DNA("ATCG")
    dna.type = DNAType.CIRCULAR
    assert dna.type == DNAType.CIRCULAR


def test_dna_name_setter() -> None:
    """
    Test that the name of a DNA object can be set correctly.
    """
    dna = DNA("ATCG")
    dna.name = "test"
    assert dna.name == "test"


def test_dna_direction_setter() -> None:
    """
    Test that the direction of a DNA sequence can be set to REVERSE.
    """
    dna = DNA("ATCG")
    dna.direction = DNADirection.REVERSE
    assert dna.direction == DNADirection.REVERSE


def test_dna_lower() -> None:
    """
    Test that the `lower` method of the `DNA` class.
    """
    dna = DNA("ATCG")
    assert dna.lower().sequence == "atcg"


def test_dna_upper() -> None:
    """
    Test the `upper` method of the `DNA` class.

    The `upper` method should return a new `DNA` object with the sequence in
    all uppercase letters.
    """
    dna = DNA("atcg")
    assert dna.upper().sequence == "ATCG"


def test_dna_complement() -> None:
    """
    Test the complement method of the DNA class.

    Creates a DNA object with the sequence "ATCG" and checks that the complement
    method returns a DNA object with the sequence "TAGC".
    """
    dna = DNA("ATCG")
    assert dna.complement().sequence == "TAGC"


def test_dna_reverse() -> None:
    """
    Test the reverse method of the DNA class.

    Creates a DNA object with the sequence "ATCG", calls the reverse method, and
    asserts that the resulting sequence is "GCTA".
    """
    dna = DNA("ATCG")
    assert dna.reverse().sequence == "GCTA"


def test_dna_eq() -> None:
    """
    Test that the equality operator works as expected for DNA sequences.
    """
    dna1 = DNA("ATCG")
    dna2 = DNA("ATCG")
    assert dna1 == dna2


def test_dna_is_complement_of() -> None:
    """
    Test that the is_complement_of method of the DNA class.
    """
    dna1 = DNA("ATCG")
    dna2 = DNA("TAGC", direction=DNADirection.REVERSE)
    assert dna1.is_complement_of(dna2)


def test_dna_len() -> None:
    """
    Test the length of a DNA sequence.

    Creates a DNA object with the sequence "ATCG" and checks that its length is 4.
    """
    dna = DNA("ATCG")
    assert len(dna) == 4


def test_dna_pad() -> None:
    """
    Test the pad method of the DNA class.

    Creates a DNA object with the sequence "ATCG" and pads it with 2 additional
    nucleotides. Asserts that the resulting sequence is "CGATCG".
    """
    dna = DNA("ATCG", dna_type=DNAType.CIRCULAR)
    assert dna.pad(2).sequence == "CGATCG"
    dna = DNA("ATCG")
    assert dna.pad(2).sequence == "--ATCG"


def test_dna_getitem() -> None:
    """
    Test the __getitem__ method of the DNA class.

    Creates a DNA object with the sequence "ATCG" and checks that slicing the
    sequence from index 1 to 3 returns a new DNA object with the sequence "TC".
    """
    dna = DNA("ATCG")
    assert dna[1:3].sequence == "TC"


def test_dna_str() -> None:
    """
    Test the string representation of a DNA object.

    Creates a DNA object with the sequence "ATCG" and checks that its string
    representation matches the expected value.
    """
    dna = DNA("ATCG")
    assert str(dna) == "DNA: ATCG"


def test_dna_repr() -> None:
    """
    Test the __repr__ method of the DNA class.
    """
    dna = DNA("ATCG")
    assert (
        repr(dna) == "DNA: Name: ATCG, DNAType.LINEAR, DNADirection.FORWARD, Seq: ATCG"
    )
