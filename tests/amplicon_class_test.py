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

"""Tests for the Amplicon class."""

import pytest

from amplifyp.amplicon import Amplicon, AmpliconGenerator
from amplifyp.dna import DNA, DNADirection, DNAType, Primer
from amplifyp.errors import (
    InvalidAmpliconRangeError,
    InvalidEndDirectionError,
    InvalidIndexOrderError,
    InvalidStartDirectionError,
    TemplateMismatchError,
)
from amplifyp.repliconf import DirIdx, Repliconf


def create_dummy_amplicon(q_score: float, circular: bool = False) -> Amplicon:
    """Helper to create a dummy Amplicon with a specific Q-score."""
    # We need valid objects for the other attributes, even if they aren't used
    # in q_score_report
    dummy_dna = DNA("A", name="dummy")
    dummy_primer = Primer("A", name="dummy_primer")

    # Valid indices for linear: start (FWD) < end (REV)
    # Valid indices for circular: start (FWD) can be > end (REV) if wrapping
    # But simpler to just use 0 -> 1 for linear to avoid validation errors

    start_idx = DirIdx(DNADirection.FWD, 0)
    end_idx = DirIdx(DNADirection.REV, 10)

    return Amplicon(
        product=dummy_dna,
        fwd_origin=dummy_primer,
        rev_origin=dummy_primer,
        start=start_idx,
        end=end_idx,
        q_score=q_score,
        circular=circular,
    )


def test_q_score_report_good() -> None:
    """Test q_score_report for 'good' scores (< 300)."""
    amplicon = create_dummy_amplicon(299.9)
    assert amplicon.q_score_report_str() == "good"
    assert amplicon.q_score_report_str(verbose=True) == "good amplification"


def test_q_score_report_okay() -> None:
    """Test q_score_report for 'okay' scores (300 <= x < 700)."""
    amplicon = create_dummy_amplicon(300.0)
    assert amplicon.q_score_report_str() == "okay"
    assert amplicon.q_score_report_str(verbose=True) == "okay amplification"

    amplicon = create_dummy_amplicon(699.9)
    assert amplicon.q_score_report_str() == "okay"


def test_q_score_report_moderate() -> None:
    """Test q_score_report for 'moderate' scores (700 <= x < 1500)."""
    amplicon = create_dummy_amplicon(700.0)
    assert amplicon.q_score_report_str() == "moderate"
    assert amplicon.q_score_report_str(verbose=True) == "moderate amplification"

    amplicon = create_dummy_amplicon(1499.9)
    assert amplicon.q_score_report_str() == "moderate"


def test_q_score_report_weak() -> None:
    """Test q_score_report for 'weak' scores (1500 <= x < 4000)."""
    amplicon = create_dummy_amplicon(1500.0)
    assert amplicon.q_score_report_str() == "weak"
    assert (
        amplicon.q_score_report_str(verbose=True)
        == "weak amplification — might be visible on an agarose gel"
    )

    amplicon = create_dummy_amplicon(3999.9)
    assert amplicon.q_score_report_str() == "weak"


def test_q_score_report_very_weak() -> None:
    """Test q_score_report for 'very weak' scores (>= 4000)."""
    amplicon = create_dummy_amplicon(4000.0)
    assert amplicon.q_score_report_str() == "very weak"
    assert (
        amplicon.q_score_report_str(verbose=True)
        == "very weak amplification — probably not visible on an agarose gel"
    )

    amplicon = create_dummy_amplicon(10000.0)
    assert amplicon.q_score_report_str() == "very weak"


def test_q_score_report_circular() -> None:
    """Test q_score_report appends '(Circular)' for circular amplicons."""
    amplicon = create_dummy_amplicon(100.0, circular=True)
    assert amplicon.q_score_report_str() == "good (Circular)"

    amplicon = create_dummy_amplicon(5000.0, circular=True)
    assert (
        amplicon.q_score_report_str(verbose=True)
        == "very weak amplification — probably not visible on an agarose gel "
        "(Circular)"
    )


def test_amplicon_post_init_validation() -> None:
    """Test the validation logic in __post_init__."""
    dummy_dna = DNA("A", name="dummy")
    dummy_primer = Primer("A", name="dummy_primer")

    # 1. Invalid start direction (REV instead of FWD)
    with pytest.raises(
        InvalidStartDirectionError, match="Start direction must be forward"
    ):
        Amplicon(
            product=dummy_dna,
            fwd_origin=dummy_primer,
            rev_origin=dummy_primer,
            start=DirIdx(DNADirection.REV, 0),
            end=DirIdx(DNADirection.REV, 10),
            q_score=100.0,
            circular=False,
        )

    # 2. Invalid end direction (FWD instead of REV)
    with pytest.raises(
        InvalidEndDirectionError, match="End direction must be reverse"
    ):
        Amplicon(
            product=dummy_dna,
            fwd_origin=dummy_primer,
            rev_origin=dummy_primer,
            start=DirIdx(DNADirection.FWD, 0),
            end=DirIdx(DNADirection.FWD, 10),
            q_score=100.0,
            circular=False,
        )

    # 3. Invalid indices for linear DNA (start > end)
    with pytest.raises(
        InvalidIndexOrderError,
        match="End index must be greater than start index for linear DNA",
    ):
        Amplicon(
            product=dummy_dna,
            fwd_origin=dummy_primer,
            rev_origin=dummy_primer,
            start=DirIdx(DNADirection.FWD, 10),
            end=DirIdx(DNADirection.REV, 0),
            q_score=100.0,
            circular=False,
        )

    # 4. Valid linear case
    Amplicon(
        product=dummy_dna,
        fwd_origin=dummy_primer,
        rev_origin=dummy_primer,
        start=DirIdx(DNADirection.FWD, 0),
        end=DirIdx(DNADirection.REV, 10),
        q_score=100.0,
        circular=False,
    )

    # 5. Valid circular case (start > end is allowed)
    Amplicon(
        product=dummy_dna,
        fwd_origin=dummy_primer,
        rev_origin=dummy_primer,
        start=DirIdx(DNADirection.FWD, 10),
        end=DirIdx(DNADirection.REV, 0),
        q_score=100.0,
        circular=True,
    )


def test_amplicon_generator_add_repliconf_errors() -> None:
    """Test error conditions when adding a Repliconf to AmpliconGenerator."""
    dna1 = DNA("AAAA", name="dna1")
    dna2 = DNA("TTTT", name="dna2")

    generator = AmpliconGenerator(dna1)
    primer = Primer("A")

    # 1. Different template
    repliconf_diff_template = Repliconf(dna2, primer)
    with pytest.raises(
        TemplateMismatchError,
        match="The Repliconf contains a different template",
    ):
        generator.add_repliconf(repliconf_diff_template)


def test_construct_amplicon_sequence_linear_invalid_order() -> None:
    """Test unreachable error in _construct_amplicon_sequence for linear DNA.

    This is hard to reach via public API because get_amplicons checks
    start < end for linear DNA. We call the private method directly.
    """
    dna = DNA("ATGC", DNAType.LINEAR)
    generator = AmpliconGenerator(dna)
    primer = Primer("A")
    repliconf = Repliconf(dna, primer)

    start_idx = DirIdx(DNADirection.FWD, 3)
    # end == start -> should trigger NotImplementedError/else block
    end_idx = DirIdx(DNADirection.REV, 3)

    with pytest.raises(
        InvalidAmpliconRangeError, match="Attempted to search for an amplicon"
    ):
        generator._construct_amplicon_sequence(
            repliconf, repliconf, start_idx, end_idx
        )


def test_construct_amplicon_sequence_linear_start_gt_end() -> None:
    """Test start > end on linear DNA returns None (pass block)."""
    dna = DNA("ATGC", DNAType.LINEAR)
    generator = AmpliconGenerator(dna)
    primer = Primer("A")
    repliconf = Repliconf(dna, primer)

    start_idx = DirIdx(DNADirection.FWD, 4)
    end_idx = DirIdx(DNADirection.REV, 3)

    seq, circular = generator._construct_amplicon_sequence(
        repliconf, repliconf, start_idx, end_idx
    )
    assert seq is None
    assert circular is False
