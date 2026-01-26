# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Tests related to settings.py."""

import pytest

from amplifyp.settings import BasePairWeightsTbl, LengthWiseWeightTbl


def test_length_wise_weight_tbl() -> None:
    """Test the LengthWiseWeightTbl class."""
    weight_tbl = LengthWiseWeightTbl(0.5, {0: 0.6, 1: 0.7})
    assert weight_tbl[0] == pytest.approx(0.6)
    assert weight_tbl[1] == pytest.approx(0.7)
    weight_tbl[2] = 0.8
    assert weight_tbl[2] == pytest.approx(0.8)
    assert weight_tbl[100] == pytest.approx(0.5)


def test_empty_length_wise_weight_tbl() -> None:
    """Test the initialisation of an empty LengthWiseWeightTbl."""
    a = LengthWiseWeightTbl()
    assert a[0] == 0


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


def test_base_pair_weights_tbl() -> None:
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
    assert npwt[("A", "C")] == pytest.approx(0.5)
    assert npwt[("G", "T")] == pytest.approx(0.5)
    assert npwt[("C", "A")] == pytest.approx(0.5)
    assert npwt[("T", "G")] == pytest.approx(0.5)
    assert npwt[("A", "A")] == pytest.approx(1.0)
    assert npwt[("C", "C")] == pytest.approx(1.0)
    assert npwt[("G", "G")] == pytest.approx(1.0)
    assert npwt[("T", "T")] == pytest.approx(1.0)

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


def test_base_pair_weights_tbl_case_sensitivity() -> None:
    """Test that BasePairWeightsTbl handles case sensitivity correctly."""
    row = "acgt-"
    col = "acgt-"
    weights = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]

    # Initialize with lowercase
    tbl = BasePairWeightsTbl(row, col, weights)

    # Check if row() and column() return uppercase (accessing row and col)
    assert tbl.row() == "ACGT"
    assert tbl.column() == "ACGT"

    # Check if __getitem__ works with uppercase keys
    assert tbl["A", "A"] == pytest.approx(1.0)

    # Check if __getitem__ works with lowercase keys
    assert tbl["a", "a"] == pytest.approx(1.0)

    # Check mixed case
    assert tbl["A", "a"] == pytest.approx(1.0)
