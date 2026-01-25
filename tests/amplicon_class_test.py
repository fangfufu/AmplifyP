"""Tests for the Amplicon class."""

import pytest

from amplifyp.amplicon import Amplicon
from amplifyp.dna import DNA, DNADirection, Primer
from amplifyp.repliconf import DirIdx


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
    with pytest.raises(ValueError, match="Start direction must be forward"):
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
    with pytest.raises(ValueError, match="End direction must be reverse"):
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
        ValueError,
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
