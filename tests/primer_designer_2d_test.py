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

"""Unit tests for the PrimerDesigner2D class and PrimerDimers2D."""

from __future__ import annotations

import pytest

from amplifyp.dimer import PrimerDimer
from amplifyp.dna import DNA
from amplifyp.primer_designer_2d import (
    FilterMetric,
    PrimerDesigner2D,
    PrimerDimers2D,
)


def test_primer_designer_2d_basic_truncation() -> None:
    """Test 2D truncation generates grid of steps with 4 dimers each."""
    fwd_dna = DNA("ATGCGTACGT")  # length 10 -> down to 8 (3 lengths: 10, 9, 8)
    rev_dna = DNA("CGTACGTACG")  # length 10 -> down to 9 (2 lengths: 10, 9)
    designer = PrimerDesigner2D(
        fwd_dna=fwd_dna,
        fwd_min_length=8,
        rev_dna=rev_dna,
        rev_min_length=9,
    )

    assert designer.fwd_dna == fwd_dna
    assert designer.fwd_min_length == 8
    assert designer.rev_dna == rev_dna
    assert designer.rev_min_length == 9

    steps = designer.all_steps
    assert len(steps) == 6  # 3 forward x 2 reverse = 6 steps

    for step in steps:
        assert isinstance(step, PrimerDimers2D)
        assert isinstance(step.fwd_fwd, PrimerDimer)
        assert isinstance(step.rev_rev, PrimerDimer)
        assert isinstance(step.fwd_rev, PrimerDimer)
        assert isinstance(step.rev_fwd, PrimerDimer)

        # Slot 1 & 2 self-dimers
        assert (
            step.fwd_fwd.primer_1.seq.upper()
            == step.fwd_fwd.primer_2.seq.upper()
        )
        assert (
            step.rev_rev.primer_1.seq.upper()
            == step.rev_rev.primer_2.seq.upper()
        )

        # Slot 3 & 4 cross-dimers
        assert (
            step.fwd_rev.primer_1.seq.upper()
            == step.fwd_fwd.primer_1.seq.upper()
        )
        assert (
            step.fwd_rev.primer_2.seq.upper()
            == step.rev_rev.primer_1.seq.upper()
        )

        assert (
            step.rev_fwd.primer_1.seq.upper()
            == step.rev_rev.primer_1.seq.upper()
        )
        assert (
            step.rev_fwd.primer_2.seq.upper()
            == step.fwd_fwd.primer_1.seq.upper()
        )


def test_primer_designer_2d_directional_truncation() -> None:
    """Test forward truncates 3' end and reverse truncates 5' end."""
    fwd_dna = DNA("ATGCGTACGT")  # length 10 -> down to 9
    rev_dna = DNA("CGTACGTACG")  # length 10 -> down to 9
    designer = PrimerDesigner2D(fwd_dna, 9, rev_dna, 9)

    steps = designer.all_steps
    # Step 0 (len 10, len 10)
    assert steps[0].fwd_fwd.primer_1.seq.upper() == "ATGCGTACGT"
    assert steps[0].rev_rev.primer_1.seq.upper() == "CGTACGTACG"

    # Step 3 (len 9 forward chopped 3' end, len 9 reverse chopped 5' end)
    assert steps[3].fwd_fwd.primer_1.seq.upper() == "ATGCGTACG"
    assert steps[3].rev_rev.primer_1.seq.upper() == "GTACGTACG"


def test_primer_designer_2d_step_metrics() -> None:
    """Test step quality and overlap metrics on PrimerDimers2D."""
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")
    designer = PrimerDesigner2D(fwd_dna, 9, rev_dna, 9)

    step = designer[0]
    dimers = [step.fwd_fwd, step.rev_rev, step.fwd_rev, step.rev_fwd]

    expected_max_q = max(d.quality for d in dimers)
    expected_mean_q = sum(d.quality for d in dimers) / 4.0
    expected_max_o = max(d.overlap for d in dimers)
    expected_mean_o = sum(d.overlap for d in dimers) / 4.0

    assert step.max_quality == expected_max_q
    assert step.mean_quality == pytest.approx(expected_mean_q)
    assert step.max_overlap == expected_max_o
    assert step.mean_overlap == pytest.approx(expected_mean_o)


def test_primer_designer_2d_invalid_inputs() -> None:
    """Test error handling for invalid sequence lengths."""
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")

    with pytest.raises(
        ValueError, match="Forward target length n must be greater than 0"
    ):
        PrimerDesigner2D(fwd_dna, 0, rev_dna, 5)

    with pytest.raises(
        ValueError, match=r"Forward target length n .* cannot exceed"
    ):
        PrimerDesigner2D(fwd_dna, 15, rev_dna, 5)

    with pytest.raises(
        ValueError, match="Reverse target length n must be greater than 0"
    ):
        PrimerDesigner2D(fwd_dna, 5, rev_dna, 0)

    with pytest.raises(
        ValueError, match=r"Reverse target length n .* cannot exceed"
    ):
        PrimerDesigner2D(fwd_dna, 5, rev_dna, 15)


def test_primer_designer_2d_filtering_max_and_mean() -> None:
    """Test atomic step filtering by MAX and MEAN metrics."""
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")

    # Filter with MAX metric
    designer_max = PrimerDesigner2D(
        fwd_dna,
        8,
        rev_dna,
        8,
        threshold=120.0,
        filter_metric=FilterMetric.MAX,
    )
    assert designer_max.filter_metric == FilterMetric.MAX
    for step in designer_max.all_steps:
        assert step.max_quality <= 120.0

    # Filter with MEAN metric
    designer_mean = PrimerDesigner2D(
        fwd_dna,
        8,
        rev_dna,
        8,
        threshold=100.0,
        filter_metric=FilterMetric.MEAN,
    )
    assert designer_mean.filter_metric == FilterMetric.MEAN
    for step in designer_mean.all_steps:
        assert step.mean_quality <= 100.0


def test_primer_designer_2d_query_methods_and_protocols() -> None:
    """Test indexing, get_step, quality_score, best_score, len, and repr."""
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")
    designer = PrimerDesigner2D(fwd_dna, 9, rev_dna, 9)

    assert len(designer) == 4
    step_0 = designer[0]
    assert designer.get_step(0) == step_0

    with pytest.raises(IndexError):
        _ = designer[10]

    # Iteration & membership
    steps_list = list(designer)
    assert len(steps_list) == 4
    assert step_0 in designer
    assert "invalid" not in designer

    # quality_score & best_score
    qs_raw = designer.quality_score(sorted=False)
    qs_sorted = designer.quality_score(sorted=True)
    assert len(qs_raw) == 4
    assert len(qs_sorted) == 4
    assert qs_sorted[0] == designer.best_score

    # String representations
    expected_repr = (
        f"PrimerDesigner2D(fwd_dna={fwd_dna!r}, fwd_min_length=9, "
        f"rev_dna={rev_dna!r}, rev_min_length=9)"
    )
    assert repr(designer) == expected_repr
    assert (
        str(designer)
        == "PrimerDesigner2D(4 steps, fwd_min_length=9, rev_min_length=9)"
    )
