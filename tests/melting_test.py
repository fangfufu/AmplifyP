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

"""Tests for the melting module."""

from dataclasses import replace

from amplifyp.melting import calculate_tm
from amplifyp.settings import DEFAULT_MELTING_SETTINGS


def test_calculate_tm_standard_sequences() -> None:
    """Test Tm calculation for standard sequences."""
    settings = DEFAULT_MELTING_SETTINGS

    # T7 Promoter: TAATACGACTCACTATAGGG (20bp)
    # Expected: ~50-52 C with 50mM Na+ and 1.5mM Mg++ (Owczarzy 2008)
    # Was ~44-45 C with Na+ only.
    t7_seq = "TAATACGACTCACTATAGGG"
    tm_t7 = calculate_tm(t7_seq, settings)
    assert 48.0 < tm_t7 < 55.0

    # M13 Fwd: CGCCAGGGTTTTCCCAGTCACGAC (24bp)
    # High GC content, expected higher Tm
    m13_seq = "CGCCAGGGTTTTCCCAGTCACGAC"
    tm_m13 = calculate_tm(m13_seq, settings)
    assert tm_m13 > tm_t7
    assert 55.0 < tm_m13 < 70.0


def test_calculate_tm_edge_cases() -> None:
    """Test Tm calculation for edge cases."""
    settings = DEFAULT_MELTING_SETTINGS

    # Empty sequence
    assert calculate_tm("", settings) == 0.0

    # Single base
    assert calculate_tm("A", settings) == 0.0

    # Two bases (minimum for NN) - effectively doesn't bind at room temp
    # Tm ~ -200C is expected for very short/unstable sequences
    assert calculate_tm("AA", settings) < 0.0


def test_salt_dependence() -> None:
    """Test that Tm increases with higher salt concentration."""
    low_salt = replace(DEFAULT_MELTING_SETTINGS, monovalent_salt_conc=50.0)
    high_salt = replace(DEFAULT_MELTING_SETTINGS, monovalent_salt_conc=1000.0)

    seq = "TAATACGACTCACTATAGGG"

    tm_low = calculate_tm(seq, low_salt)
    tm_high = calculate_tm(seq, high_salt)

    # Higher salt should stabilize DNA, increasing Tm
    assert tm_high > tm_low


def test_magnesium_stabilization() -> None:
    """Test that adding Mg2+ increases Tm."""
    # 50 mM Na+, 0 mM Mg++
    no_mg = replace(DEFAULT_MELTING_SETTINGS, divalent_salt_conc=0.0)
    # 50 mM Na+, 1.5 mM Mg++
    with_mg = replace(DEFAULT_MELTING_SETTINGS, divalent_salt_conc=1.5)

    seq = "TAATACGACTCACTATAGGG"

    tm_no_mg = calculate_tm(seq, no_mg)
    tm_with_mg = calculate_tm(seq, with_mg)

    assert tm_with_mg > tm_no_mg
    # Expect roughly 5-10C increase
    assert (tm_with_mg - tm_no_mg) > 2.0


def test_gc_content_contribution() -> None:
    """Test that GC-rich sequences have higher Tm than AT-rich sequences.

    Checks sequences of same length.
    """
    settings = DEFAULT_MELTING_SETTINGS

    poly_a = "AAAAAAAAAAAAAAAAAAAA"  # 20 A
    poly_c = "CCCCCCCCCCCCCCCCCCCC"  # 20 C

    tm_a = calculate_tm(poly_a, settings)
    tm_c = calculate_tm(poly_c, settings)

    assert tm_c > tm_a


def test_invalid_chars_handling() -> None:
    """Test handling of invalid characters.

    Can handle skipped/ignored characters without error.
    """
    settings = DEFAULT_MELTING_SETTINGS
    # NN should only skip the invalid dinucleotide steps
    # "ACGT" -> AC, CG, GT
    # "ACNRT" -> AC, CN(skip), NR(skip), RT(skip)?
    # Current implementation: if dinuc not in table, pass.
    # So "ACNT" -> AC (valid), CN (invalid), NT(invalid). Only first step
    # counts.

    seq_normal = "ACGT"
    seq_n = "ACNT"

    tm_normal = calculate_tm(seq_normal, settings)
    tm_n = calculate_tm(seq_n, settings)

    # With N, we lose detailed energy but shouldn't crash.
    # The resulting Tm will be very low due to missing stacking energy.
    assert tm_n < tm_normal
