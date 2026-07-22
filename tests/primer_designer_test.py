# Copyright (C) 2026 AmplifyP Contributors
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

"""Unit tests for the PrimerDesigner1D class."""

from __future__ import annotations

import pytest

from amplifyp.dimer import PrimerDimer
from amplifyp.dna import DNA, DNADirection
from amplifyp.primer_designer import PrimerDesigner1D


def test_primer_designer_forward_mode() -> None:
    """Test Forward mode truncates from 3' end down to min_length."""
    dna_obj = DNA("ATGCGTACGT")  # length 10
    min_length = 7
    designer = PrimerDesigner1D(dna_obj, min_length, mode=DNADirection.FWD)

    assert designer.dna == dna_obj
    assert designer.min_length == 7
    assert designer.mode == DNADirection.FWD

    results = designer.results
    assert len(results) == 4  # lengths 10, 9, 8, 7

    assert [len(dimer.primer_1) for dimer in results] == [10, 9, 8, 7]
    assert [dimer.primer_1.seq.upper() for dimer in results] == [
        "ATGCGTACGT",
        "ATGCGTACG",
        "ATGCGTAC",
        "ATGCGTA",
    ]

    for dimer in results:
        assert isinstance(dimer, PrimerDimer)
        assert dimer.primer_1.seq.upper() == dimer.primer_2.seq.upper()


def test_primer_designer_backward_mode() -> None:
    """Test Backward mode truncates from 5' end down to min_length."""
    dna_obj = DNA("ATGCGTACGT")  # length 10
    min_length = 7
    designer = PrimerDesigner1D(dna_obj, min_length, mode=DNADirection.REV)

    assert designer.dna == dna_obj
    assert designer.min_length == 7
    assert designer.mode == DNADirection.REV

    results = designer.results
    assert len(results) == 4  # lengths 10, 9, 8, 7

    assert [len(dimer.primer_1) for dimer in results] == [10, 9, 8, 7]
    assert [dimer.primer_1.seq.upper() for dimer in results] == [
        "ATGCGTACGT",
        "TGCGTACGT",
        "GCGTACGT",
        "CGTACGT",
    ]


def test_primer_designer_default_mode() -> None:
    """Test default mode is DNADirection.FWD."""
    dna_obj = DNA("ATGCGTACGT")
    designer = PrimerDesigner1D(dna_obj, 8)
    assert designer.mode == DNADirection.FWD


def test_primer_designer_query_methods() -> None:
    """Test indexing, len, get_dimer, and best_score query methods."""
    dna_obj = DNA(
        "ATGCGTACGT"
    )  # length 10 down to 7 -> 4 steps (index 0, 1, 2, 3)
    designer = PrimerDesigner1D(dna_obj, 7, mode=DNADirection.FWD)

    assert len(designer) == 4

    dimer_0 = designer[0]
    assert len(dimer_0.primer_1) == 10
    assert dimer_0.primer_1.seq.upper() == "ATGCGTACGT"

    dimer_2 = designer[2]  # length 8 step
    assert len(dimer_2.primer_1) == 8
    assert dimer_2.primer_1.seq.upper() == "ATGCGTAC"

    assert designer.get_dimer(2) == dimer_2

    with pytest.raises(IndexError):
        _ = designer[10]

    best_pair = designer.best_score
    assert isinstance(best_pair, tuple)
    assert len(best_pair) == 2
    best_idx, best_q = best_pair
    assert designer[best_idx] == min(designer.results, key=lambda d: d.quality)
    assert best_q == designer[best_idx].quality


def test_primer_designer_invalid_inputs() -> None:
    """Test error handling for invalid sequence lengths, modes, and types."""
    dna_obj = DNA("ATGCGTACGT")

    with pytest.raises(
        ValueError, match="Target length n must be greater than 0"
    ):
        PrimerDesigner1D(dna_obj, 0)

    with pytest.raises(
        ValueError, match="cannot exceed initial sequence length"
    ):
        PrimerDesigner1D(dna_obj, 15)


def test_primer_designer_sequence_protocol_and_representations() -> None:
    """Test iteration, membership, repr, str, and query additions."""
    dna_obj = DNA("ATGCGTACGT")
    designer = PrimerDesigner1D(dna_obj, min_length=7)

    assert designer.min_length == 7
    assert len(designer) == 4

    # Iteration
    dimers_list = list(designer)
    assert len(dimers_list) == 4
    assert dimers_list == list(designer.results)

    # Membership
    assert designer[0] in designer
    assert "not_a_dimer" not in designer

    # repr and str
    expected_repr = (
        f"PrimerDesigner1D(dna={dna_obj!r}, min_length=7, "
        f"mode={DNADirection.FWD!r})"
    )
    assert repr(designer) == expected_repr
    assert str(designer) == "PrimerDesigner1D(4 steps, min_length=7, mode=FWD)"

    # quality_score
    sorted_qs = designer.quality_score(sorted=True)

    assert len(sorted_qs) == 4
    assert [score for _, score in sorted_qs] == sorted(
        d.quality for d in designer.results
    )

    raw_qs = designer.quality_score(sorted=False)
    assert raw_qs == tuple(
        (i, d.quality) for i, d in enumerate(designer.results)
    )

    assert sorted_qs[0] == designer.best_score


def test_primer_designer_custom_generator() -> None:
    """Test passing custom PrimerDimerGenerator."""
    from amplifyp.dimer import PrimerDimerGenerator
    from amplifyp.settings import PrimerDimerSettings

    custom_settings = PrimerDimerSettings(min_overlap=3)
    custom_gen = PrimerDimerGenerator(settings=custom_settings)

    dna_obj = DNA("ATGCGTACGT")
    designer = PrimerDesigner1D(dna_obj, 7, generator=custom_gen)

    assert designer.generator == custom_gen
    assert designer[0].settings == custom_settings
