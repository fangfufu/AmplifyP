# -*- coding: utf-8 -*-
# pylint: disable=cyclic-import
"""Tests for types.py."""

import pytest

from amplifyp.dna import (
    DNA,
    DNAType,
    Nucleotides,
)

from amplifyp.settings import (
    LengthWiseWeightTbl,
    BasePairWeights,
)


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


def test_run_length_weight_tbl() -> None:
    """Test the LengthWiseWeightTbl class."""
    weight_tbl = LengthWiseWeightTbl(5, 0.5, [(0, 0.6), (1, 0.7)])
    assert weight_tbl[0] == 0.6
    assert weight_tbl[1] == 0.7
    weight_tbl[2] = 0.8
    assert weight_tbl[2] == 0.8
    assert list(weight_tbl) == [0.6, 0.7, 0.8, 0.5, 0.5]


def test_nucleotide_pairwise_weight_tbl() -> None:
    """Test function for the NucleotidePairwiseWeightTbl class."""
    # Define the nucleotides and their pairwise weights
    nucleotides = "ACGT"
    pairwise_weights = [
        [1.0, 0.5, 0.2, 0.1],
        [0.5, 1.0, 0.1, 0.2],
        [0.2, 0.1, 0.2, 0.5],
        [0.1, 0.2, 0.5, 0.2],
    ]

    # Create a NucleotidePairwiseWeightTbl instance
    npwt = BasePairWeights(nucleotides, nucleotides, pairwise_weights)

    # Test the row and column properties
    assert npwt.row == ["A", "C", "G", "T"]
    assert npwt.column == ["A", "C", "G", "T"]

    # Test the __setitem__ method
    npwt[("G", "G")] = 1.0
    npwt[("T", "T")] = 1.0

    # Test the __getitem__ method
    assert npwt[("A", "C")] == 0.5
    assert npwt[("G", "T")] == 0.5
    assert npwt[("C", "A")] == 0.5
    assert npwt[("T", "G")] == 0.5
    assert npwt[("A", "A")] == 1.0
    assert npwt[("C", "C")] == 1.0
    assert npwt[("G", "G")] == 1.0
    assert npwt[("T", "T")] == 1.0

    # Create an invalid NucleotidePairwiseWeightTbl instance
    with pytest.raises(ValueError):
        npwt = BasePairWeights(Nucleotides.PRIMER, Nucleotides.TARGET, pairwise_weights)
