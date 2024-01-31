# -*- coding: utf-8 -*-
"""Tests related to settings.py."""

import pytest
from amplifyp.settings import BasePairWeightsTbl, LengthWiseWeightTbl, DEFAULT_SETTINGS


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


pairwise_weights = [
    [1.0, 0.5, 0.2, 0.1],
    [0.5, 1.0, 0.1, 0.2],
    [0.2, 0.1, 0.2, 0.5],
    [0.1, 0.2, 0.5, 0.2],
]


def test_invalid_tbl_generation() -> None:
    """Test invalid BasePairWeightsTbl creation."""

    with pytest.raises(ValueError):
        BasePairWeightsTbl("AB", "ABCD", pairwise_weights)

    with pytest.raises(ValueError):
        BasePairWeightsTbl("ABCD", "AB", pairwise_weights)


def test_nucleotide_pairwise_weight_tbl() -> None:
    """Test function for the BasePairWeightsTbl class."""
    # Define the nucleotides and their pairwise weights
    nucleotides = "ACGT-"

    # Create a BasePairWeightsTbl instance
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

    # Test the support for gap symbol
    assert npwt["-", "A"] == 0

    assert (
        str(npwt)
        == "{('A', 'A'): 1.0, ('A', 'C'): 0.5, ('A', 'G'): 0.2, \
('A', 'T'): 0.1, ('A', '-'): 0, ('C', 'A'): 0.5, ('C', 'C'): 1.0, ('C', 'G'): \
0.1, ('C', 'T'): 0.2, ('C', '-'): 0, ('G', 'A'): 0.2, ('G', 'C'): 0.1, \
('G', 'G'): 1.0, ('G', 'T'): 0.5, ('G', '-'): 0, ('T', 'A'): 0.1, ('T', 'C'): \
0.2, ('T', 'G'): 0.5, ('T', 'T'): 1.0, ('T', '-'): 0, ('-', 'A'): 0, \
('-', 'C'): 0, ('-', 'G'): 0, ('-', 'T'): 0, ('-', '-'): 0}"
    )

    assert len(npwt) == 16
