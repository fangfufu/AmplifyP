# -*- coding: utf-8 -*-
"""Tests for types.py."""
from amplifyp.types import NucleotidePairwiseWeightTbl, RunLengthWeightTbl


def test_run_length_weight_tbl() -> None:
    """Test the RunLengthWeightTbl class."""
    weight_tbl = RunLengthWeightTbl(5, 0.5)
    weight_tbl[0] = 0.6
    weight_tbl[1] = 0.7
    weight_tbl[2] = 0.8
    weight_tbl[3] = 0.9
    weight_tbl[4] = 1.0
    assert str(weight_tbl) == "[0.6, 0.7, 0.8, 0.9, 1.0]"
    assert repr(weight_tbl) == "[0.6, 0.7, 0.8, 0.9, 1.0]"
    assert list(weight_tbl) == [0.6, 0.7, 0.8, 0.9, 1.0]


def test_nucleotide_pairwise_weight_tbl() -> None:
    """Test function for the NucleotidePairwiseWeightTbl class.

    This function tests the initialization of the
    NucleotidePairwiseWeightTbl class, setting and getting values, and
    iteration over the table.
    """
    # Test initialization
    tbl = NucleotidePairwiseWeightTbl("ACGT", "ACGT", 1.0)
    assert len(tbl) == 16
    assert tbl[("A", "A")] == 1.0
    assert tbl[("A", "C")] == 1.0
    assert tbl[("C", "A")] == 1.0
    assert tbl[("C", "C")] == 1.0

    # Test setting and getting values
    tbl[("A", "A")] = 2.0
    assert tbl[("A", "A")] == 2.0

    # Test iteration
    keys = set()
    for key in tbl:
        keys.add(key)
    assert keys == {
        ("A", "A"),
        ("A", "C"),
        ("A", "G"),
        ("A", "T"),
        ("C", "A"),
        ("C", "C"),
        ("C", "G"),
        ("C", "T"),
        ("G", "A"),
        ("G", "C"),
        ("G", "G"),
        ("G", "T"),
        ("T", "A"),
        ("T", "C"),
        ("T", "G"),
        ("T", "T"),
    }
