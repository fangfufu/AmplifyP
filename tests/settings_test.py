# -*- coding: utf-8 -*-
"""Tests for types.py."""

import pytest
from amplifyp.settings import BasePairWeightsTbl, LengthWiseWeightTbl, DEFAULT_SETTINGS
from amplifyp.dna import Nucleotides


def test_run_length_weight_tbl() -> None:
    """Test the LengthWiseWeightTbl class."""
    weight_tbl = LengthWiseWeightTbl(0.5, {0: 0.6, 1: 0.7})
    assert weight_tbl[0] == 0.6
    assert weight_tbl[1] == 0.7
    weight_tbl[2] = 0.8
    assert weight_tbl[2] == 0.8
    assert weight_tbl[100] == 0.5


def test_default_run_weight_tbl() -> None:
    """Test if the default settings work"""
    assert DEFAULT_SETTINGS.match_weight[0] == 30


def test_nucleotide_pairwise_weight_tbl() -> None:
    """Test function for the NucleotidePairwiseWeightTbl class."""
    # Define the nucleotides and their pairwise weights
    nucleotides = "ACGT-"
    pairwise_weights = [
        [1.0, 0.5, 0.2, 0.1],
        [0.5, 1.0, 0.1, 0.2],
        [0.2, 0.1, 0.2, 0.5],
        [0.1, 0.2, 0.5, 0.2],
    ]

    # Create a NucleotidePairwiseWeightTbl instance
    npwt = BasePairWeightsTbl(nucleotides, nucleotides, pairwise_weights)

    # Test the row and column properties
    assert npwt.row() == "ACGT"
    assert npwt.column() == "ACGT"

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
        npwt = BasePairWeightsTbl(
            Nucleotides.PRIMER, Nucleotides.LINEAR, pairwise_weights
        )

    # Test the support for gap symbol
    assert npwt["-", "A"] == 0
