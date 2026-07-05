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

"""Tests for the melting module."""

from dataclasses import replace

import pytest

from amplifyp.dna import Primer
from amplifyp.errors import InsufficientThermodynamicDataError
from amplifyp.melting import (
    calculate_tm_lander_amplify4,
    calculate_tm_santalucia_1998_owczarzy_2008,
)
from amplifyp.settings import GLOBAL_TM_SETTINGS
from tests.examples.amplify4_examples import (
    primer_11bp,
    primer_1701,
    primer_10289,
    primer_10290,
)


def test_calculate_tm_standard_sequences() -> None:
    """Test Tm calculation for standard sequences."""
    settings = GLOBAL_TM_SETTINGS

    # T7 Promoter: TAATACGACTCACTATAGGG (20bp)
    # Expected: ~50-52 C with 50mM Na+ and 1.5mM Mg++ (Owczarzy 2008)
    # Was ~44-45 C with Na+ only.
    t7_seq = "TAATACGACTCACTATAGGG"
    tm_t7 = calculate_tm_santalucia_1998_owczarzy_2008(Primer(t7_seq), settings)
    assert 48.0 < tm_t7 < 55.0

    # M13 Fwd: CGCCAGGGTTTTCCCAGTCACGAC (24bp)
    # High GC content, expected higher Tm
    m13_seq = "CGCCAGGGTTTTCCCAGTCACGAC"
    tm_m13 = calculate_tm_santalucia_1998_owczarzy_2008(
        Primer(m13_seq), settings
    )
    assert tm_m13 > tm_t7
    assert 55.0 < tm_m13 < 70.0


def test_calculate_tm_edge_cases() -> None:
    """Test Tm calculation for edge cases."""
    settings = GLOBAL_TM_SETTINGS

    # Empty sequence
    assert calculate_tm_santalucia_1998_owczarzy_2008(
        Primer(""), settings
    ) == pytest.approx(0.0)

    # Single base
    assert calculate_tm_santalucia_1998_owczarzy_2008(
        Primer("A"), settings
    ) == pytest.approx(0.0)

    # Two bases (minimum for NN) - effectively doesn't bind at room temp
    # Tm ~ -200C is expected for very short/unstable sequences
    assert (
        calculate_tm_santalucia_1998_owczarzy_2008(Primer("AA"), settings) < 0.0
    )


def test_salt_dependence() -> None:
    """Test that Tm increases with higher salt concentration."""
    low_salt = replace(GLOBAL_TM_SETTINGS, monovalent_salt_conc=50.0)
    high_salt = replace(GLOBAL_TM_SETTINGS, monovalent_salt_conc=1000.0)

    seq = "TAATACGACTCACTATAGGG"

    tm_low = calculate_tm_santalucia_1998_owczarzy_2008(Primer(seq), low_salt)
    tm_high = calculate_tm_santalucia_1998_owczarzy_2008(Primer(seq), high_salt)

    # Higher salt should stabilise DNA, increasing Tm
    assert tm_high > tm_low


def test_magnesium_stabilisation() -> None:
    """Test that adding Mg2+ increases Tm."""
    # 50 mM Na+, 0 mM Mg++
    no_mg = replace(GLOBAL_TM_SETTINGS, divalent_salt_conc=0.0)
    # 50 mM Na+, 1.5 mM Mg++
    with_mg = replace(GLOBAL_TM_SETTINGS, divalent_salt_conc=1.5)

    seq = "TAATACGACTCACTATAGGG"

    tm_no_mg = calculate_tm_santalucia_1998_owczarzy_2008(Primer(seq), no_mg)
    tm_with_mg = calculate_tm_santalucia_1998_owczarzy_2008(
        Primer(seq), with_mg
    )

    assert tm_with_mg > tm_no_mg
    # Expect roughly 5-10C increase
    assert (tm_with_mg - tm_no_mg) > 2.0


def test_gc_content_contribution() -> None:
    """Test that GC-rich sequences have higher Tm than AT-rich sequences.

    Checks sequences of same length.
    """
    settings = GLOBAL_TM_SETTINGS

    poly_a = "AAAAAAAAAAAAAAAAAAAA"  # 20 A
    poly_c = "CCCCCCCCCCCCCCCCCCCC"  # 20 C

    tm_a = calculate_tm_santalucia_1998_owczarzy_2008(Primer(poly_a), settings)
    tm_c = calculate_tm_santalucia_1998_owczarzy_2008(Primer(poly_c), settings)

    assert tm_c > tm_a


def test_invalid_chars_handling() -> None:
    """Test handling of invalid characters.

    Can handle skipped/ignored characters without error.
    """
    settings = GLOBAL_TM_SETTINGS
    # NN should only skip the invalid dinucleotide steps
    # "ACGT" -> AC, CG, GT
    # "ACNRT" -> AC, CN(skip), NR(skip), RT(skip)?
    # Current implementation: if dinuc not in table, pass.
    # So "ACNT" -> AC (valid), CN (invalid), NT(invalid). Only first step
    # counts.

    seq_normal = "ACGT"
    seq_n = "ACNT"

    tm_normal = calculate_tm_santalucia_1998_owczarzy_2008(
        Primer(seq_normal), settings
    )
    tm_n = calculate_tm_santalucia_1998_owczarzy_2008(Primer(seq_n), settings)

    # With N, we lose detailed energy but shouldn't crash.
    # The resulting Tm will be very low due to missing stacking energy.
    assert tm_n < tm_normal


def test_amplify4_tm_calculation() -> None:
    """Test Tm calculation using Amplify4 algorithm."""
    assert calculate_tm_lander_amplify4(
        primer_11bp, GLOBAL_TM_SETTINGS
    ) == pytest.approx(19.48, abs=0.01)
    assert calculate_tm_lander_amplify4(
        primer_1701, GLOBAL_TM_SETTINGS
    ) == pytest.approx(64.01, abs=0.01)
    assert calculate_tm_lander_amplify4(
        primer_10289, GLOBAL_TM_SETTINGS
    ) == pytest.approx(69.04, abs=0.01)
    assert calculate_tm_lander_amplify4(
        primer_10290, GLOBAL_TM_SETTINGS
    ) == pytest.approx(57.64, abs=0.01)


def test_tm_calculation() -> None:
    """Test Tm calculation using standard algorithm."""
    assert calculate_tm_santalucia_1998_owczarzy_2008(
        primer_11bp, GLOBAL_TM_SETTINGS
    ) == pytest.approx(23.251668245623136)
    assert calculate_tm_santalucia_1998_owczarzy_2008(
        primer_1701, GLOBAL_TM_SETTINGS
    ) == pytest.approx(60.15678755876894)
    assert calculate_tm_santalucia_1998_owczarzy_2008(
        primer_10289, GLOBAL_TM_SETTINGS
    ) == pytest.approx(65.15445294119144)
    assert calculate_tm_santalucia_1998_owczarzy_2008(
        primer_10290, GLOBAL_TM_SETTINGS
    ) == pytest.approx(55.377546798740696)


def test_calculate_tm_zero_conc() -> None:
    """Test Tm calculation with zero DNA concentration."""
    # Standard settings have 50nM.
    # We want to test the branch where dna_conc <= 0 -> 50e-9

    zero_conc_settings = replace(GLOBAL_TM_SETTINGS, dna_conc=0.0)

    # Should give same result as 50nM if it falls back to 50nM
    # Or at least shouldn't crash
    tm = calculate_tm_santalucia_1998_owczarzy_2008(
        Primer("ATCG"), zero_conc_settings
    )
    assert tm > -273.15


def test_calculate_tm_no_monovalent() -> None:
    """Test Tm calculation with no monovalent salt."""
    # This hits the divergence check mono_M == 0
    settings = replace(GLOBAL_TM_SETTINGS, monovalent_salt_conc=0.0)
    tm = calculate_tm_santalucia_1998_owczarzy_2008(Primer("ATCG"), settings)
    assert tm > -273.15


def test_calculate_tm_amplify4_edge_cases() -> None:
    """Test Tm calculation (Amplify4) with edge cases."""
    # Empty primer
    assert calculate_tm_lander_amplify4(
        Primer(""), GLOBAL_TM_SETTINGS
    ) == pytest.approx(0.0)

    # Zero/negative salt
    # Should fallback to 50mM
    neg_salt_settings = replace(GLOBAL_TM_SETTINGS, monovalent_salt_conc=-10.0)
    tm = calculate_tm_lander_amplify4(primer_11bp, neg_salt_settings)
    assert tm > 0.0


def test_calculate_tm_no_salts() -> None:
    """Test Tm calculation with no salts (monovalent=0, divalent=0)."""
    settings = replace(
        GLOBAL_TM_SETTINGS, monovalent_salt_conc=0.0, divalent_salt_conc=0.0
    )
    # Should fallback to 1M Na+ logic (tm_1m_K) without correction?
    # Or just return tm_1m_K essentially.
    tm = calculate_tm_santalucia_1998_owczarzy_2008(Primer("ATCG"), settings)
    assert tm > -273.15


def test_calculate_tm_invalid_sequences() -> None:
    """Test that Tm calculation raises ValueError for invalid sequences."""
    settings = GLOBAL_TM_SETTINGS
    with pytest.raises(InsufficientThermodynamicDataError):
        calculate_tm_santalucia_1998_owczarzy_2008(Primer("NNNN"), settings)
