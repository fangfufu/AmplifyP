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

import threading

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


def test_primer_designer_2d_amplicon_filtering() -> None:
    """Test template DNA input and max_amplicon_count filtering."""
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("ACGTACGCAT")
    template_dna = DNA("ATGCGTACGTTTTATGCGTACGTTTTATGCGTACGT")

    designer = PrimerDesigner2D(
        fwd_dna,
        8,
        rev_dna,
        8,
        template=template_dna,
        max_amplicon_count=1,
    )
    assert designer.template == template_dna
    assert designer.max_amplicon_count == 1

    from amplifyp.amplicon import AmpliconGenerator
    from amplifyp.repliconf import Repliconf

    for step in designer.all_steps:
        amp_gen = AmpliconGenerator(template_dna)
        fwd_conf = Repliconf(template_dna, step.fwd_fwd.primer_1)
        amp_gen.add_repliconf(fwd_conf)
        if (
            step.fwd_fwd.primer_1.seq.upper()
            != step.rev_rev.primer_1.seq.upper()
        ):
            rev_conf = Repliconf(template_dna, step.rev_rev.primer_1)
            amp_gen.add_repliconf(rev_conf)
        assert len(amp_gen.get_amplicons()) <= 1
        assert step.amplicon_count is not None
        assert step.amplicon_count <= 1


def test_primer_designer_2d_template_unconstrained_amplicons() -> None:
    """Test template DNA input with unconstrained amplicons."""
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("ACGTACGCAT")
    template_dna = DNA("ATGCGTACGTTTTATGCGTACGTTTTATGCGTACGT")

    designer = PrimerDesigner2D(
        fwd_dna,
        8,
        rev_dna,
        8,
        template=template_dna,
        max_amplicon_count=None,
    )
    assert designer.template == template_dna
    assert designer.max_amplicon_count is None
    assert len(designer.all_steps) == 9  # 3 fwd x 3 rev
    for step in designer.all_steps:
        assert step.amplicon_count is not None
        assert isinstance(step.amplicon_count, int)
        assert step.amplicon_count >= 1


def test_primer_designer_2d_zero_amplicons_filtered() -> None:
    """Test that candidate pairs generating 0 amplicons are excluded."""
    fwd_dna = DNA("AAAAAAAAAA")
    rev_dna = DNA("CCCCCCCCCC")
    # Template has no binding sites for these primers
    template_dna = DNA("GGGGGGGGGGGGGGGGGGGG")

    designer = PrimerDesigner2D(
        fwd_dna,
        8,
        rev_dna,
        8,
        template=template_dna,
    )
    assert len(designer.all_steps) == 0


def test_primer_designer_2d_amplicon_invalid_inputs() -> None:
    """Test validation errors for template and max_amplicon_count."""
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")
    template_dna = DNA("ATGCGTACGTTTTATGCGTACGT")

    with pytest.raises(ValueError, match="Template DNA must be provided"):
        PrimerDesigner2D(fwd_dna, 8, rev_dna, 8, max_amplicon_count=1)

    with pytest.raises(
        ValueError, match="max_amplicon_count must be non-negative"
    ):
        PrimerDesigner2D(
            fwd_dna,
            8,
            rev_dna,
            8,
            template=template_dna,
            max_amplicon_count=-1,
        )


def test_primer_designer_2d_empty_and_mean_score() -> None:
    """Test 2D designer empty best_score error and FilterMetric.MEAN."""

    from amplifyp.dimer import PrimerDimerGenerator

    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")

    # Empty steps raises RuntimeError on best_score
    designer_empty = PrimerDesigner2D(fwd_dna, 9, rev_dna, 9, threshold=-999.0)
    with pytest.raises(RuntimeError, match="No analysis steps recorded"):
        _ = designer_empty.best_score

    # FilterMetric.MEAN
    custom_gen = PrimerDimerGenerator()
    designer_mean = PrimerDesigner2D(
        fwd_dna,
        9,
        rev_dna,
        9,
        generator=custom_gen,
        threshold=200.0,
        max_overlap=6,
        filter_metric=FilterMetric.MEAN,
    )
    assert designer_mean.generator == custom_gen
    assert designer_mean.threshold == 200.0
    assert designer_mean.max_overlap == 6
    scores = designer_mean.quality_score(sorted=True)
    assert len(scores) > 0


def test_primer_designer_2d_on_progress_callback() -> None:
    """Test on_progress callback fires once per combination with
    correct args."""
    # 3 fwd lengths (10, 9, 8) x 2 rev lengths (10, 9) = 6 combinations
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")
    fwd_min, rev_min = 8, 9
    expected_total = (len(fwd_dna.seq) - fwd_min + 1) * (
        len(rev_dna.seq) - rev_min + 1
    )

    calls: list[tuple[int, int]] = []

    def _cb(done: int, total: int) -> None:
        calls.append((done, total))

    PrimerDesigner2D(
        fwd_dna=fwd_dna,
        fwd_min_length=fwd_min,
        rev_dna=rev_dna,
        rev_min_length=rev_min,
        on_progress=_cb,
    )

    # Exactly one call per combination.
    assert len(calls) == expected_total
    # total is consistent across all calls.
    assert all(t == expected_total for _, t in calls)
    # done values are strictly increasing from 1 to total.
    done_values = [d for d, _ in calls]
    assert done_values == list(range(1, expected_total + 1))


def test_primer_designer_2d_no_cancel_event_not_aborted() -> None:
    """Test aborted is False when no cancel event is provided."""
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")
    designer = PrimerDesigner2D(
        fwd_dna=fwd_dna,
        fwd_min_length=8,
        rev_dna=rev_dna,
        rev_min_length=9,
    )
    assert designer.aborted is False


def test_primer_designer_2d_cancel_event_aborts_early() -> None:
    """Test cancel event stops analysis and retains partial results."""
    # 3 fwd lengths (10, 9, 8) x 2 rev lengths (10, 9) = 6 combinations
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")
    fwd_min, rev_min = 8, 9
    expected_total = (len(fwd_dna.seq) - fwd_min + 1) * (
        len(rev_dna.seq) - rev_min + 1
    )
    cancel = threading.Event()

    def _cb(done: int, total: int) -> None:
        if done == 3:
            cancel.set()

    designer = PrimerDesigner2D(
        fwd_dna=fwd_dna,
        fwd_min_length=fwd_min,
        rev_dna=rev_dna,
        rev_min_length=rev_min,
        on_progress=_cb,
        cancel_event=cancel,
    )

    # Analysis stopped after the combination that triggered the cancel.
    assert designer.aborted is True
    assert len(designer) == 3
    assert len(designer) < expected_total


def test_primer_designer_2d_cancel_event_pre_set() -> None:
    """Test a pre-set cancel event aborts before any step completes."""
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")
    cancel = threading.Event()
    cancel.set()

    designer = PrimerDesigner2D(
        fwd_dna=fwd_dna,
        fwd_min_length=8,
        rev_dna=rev_dna,
        rev_min_length=9,
        cancel_event=cancel,
    )

    assert designer.aborted is True
    assert len(designer) == 0


def test_primer_designer_2d_unset_cancel_event_runs_to_completion() -> None:
    """Test an unset cancel event does not affect full analysis."""
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")
    fwd_min, rev_min = 8, 9
    expected_total = (len(fwd_dna.seq) - fwd_min + 1) * (
        len(rev_dna.seq) - rev_min + 1
    )
    cancel = threading.Event()

    designer = PrimerDesigner2D(
        fwd_dna=fwd_dna,
        fwd_min_length=fwd_min,
        rev_dna=rev_dna,
        rev_min_length=rev_min,
        cancel_event=cancel,
    )

    assert designer.aborted is False
    assert len(designer) == expected_total


def test_primer_designer_2d_on_progress_skipped_branches() -> None:
    """Test progress callback invocations when steps are filtered out."""
    fwd_dna = DNA("ATGCGTACGT")
    rev_dna = DNA("CGTACGTACG")
    progress_calls: list[tuple[int, int]] = []

    def on_prog(done: int, total: int) -> None:
        progress_calls.append((done, total))

    # 1. amplicon_count < 1 with template that produces 0 amplicons
    template = DNA("A" * 50)
    designer_tpl = PrimerDesigner2D(
        fwd_dna=fwd_dna,
        fwd_min_length=10,
        rev_dna=rev_dna,
        rev_min_length=10,
        template=template,
        on_progress=on_prog,
    )
    assert len(designer_tpl) == 0
    assert len(progress_calls) == 1

    # 2. threshold filter skip with progress callback
    progress_calls.clear()
    designer_thresh = PrimerDesigner2D(
        fwd_dna=fwd_dna,
        fwd_min_length=10,
        rev_dna=rev_dna,
        rev_min_length=10,
        threshold=-1.0,
        on_progress=on_prog,
    )
    assert len(designer_thresh) == 0
    assert len(progress_calls) == 1

    # 3. max_overlap filter skip with progress callback
    progress_calls.clear()
    designer_overlap = PrimerDesigner2D(
        fwd_dna=fwd_dna,
        fwd_min_length=10,
        rev_dna=rev_dna,
        rev_min_length=10,
        max_overlap=-1,
        on_progress=on_prog,
    )
    assert len(designer_overlap) == 0
    assert len(progress_calls) == 1

    # 4. max_amplicon_count filter skip with progress callback
    progress_calls.clear()
    designer_max_amp = PrimerDesigner2D(
        fwd_dna=DNA("ATGCATGCATGC"),
        fwd_min_length=12,
        rev_dna=DNA("GCATGCATGCAT"),
        rev_min_length=12,
        template=DNA("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"),
        max_amplicon_count=0,
        on_progress=on_prog,
    )
    assert len(designer_max_amp) == 0
    assert len(progress_calls) == 1
