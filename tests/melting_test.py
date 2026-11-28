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
from typing import Any
from unittest.mock import MagicMock, patch

import pytest

from amplifyp.dna import Primer
from amplifyp.errors import InsufficientThermodynamicDataError
from amplifyp.melting import (
    IUPAC_BASES,
    TmMethod,
    TmRange,
    calculate_tm,
    calculate_tm_lander_amplify4,
    calculate_tm_range,
    calculate_tm_santalucia_1998_owczarzy_2008,
)
from amplifyp.settings import GLOBAL_TM_SETTINGS, TMSettings
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
    """Test handling of IUPAC wildcard bases and non-IUPAC rejection."""
    settings = GLOBAL_TM_SETTINGS
    seq_n = "ACNT"
    tm_n = calculate_tm_santalucia_1998_owczarzy_2008(Primer(seq_n), settings)
    assert tm_n == pytest.approx(calculate_tm(Primer(seq_n), settings))

    bad_primer = MagicMock(spec=Primer)
    bad_primer.seq_upper = "AC-T"
    with pytest.raises(InsufficientThermodynamicDataError):
        _ = calculate_tm_santalucia_1998_owczarzy_2008(bad_primer, settings)


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
    """Test that Tm calculation raises for non-IUPAC sequences and handles N."""
    settings = GLOBAL_TM_SETTINGS
    assert (
        calculate_tm_santalucia_1998_owczarzy_2008(Primer("NNNN"), settings)
        > -273.15
    )

    bad_primer = MagicMock(spec=Primer)
    bad_primer.seq_upper = "XXXX"
    with pytest.raises(InsufficientThermodynamicDataError):
        calculate_tm_santalucia_1998_owczarzy_2008(bad_primer, settings)


def test_dntp_chelation() -> None:
    """Test that dNTPs chelate Mg2+ and reduce free Mg concentration."""
    # 50 mM Na+, 1.5 mM Mg2+, 0.0 mM dNTP
    no_dntp = TMSettings(divalent_salt_conc=1.5, dntp_conc=0.0)
    # 50 mM Na+, 1.5 mM Mg2+, 0.8 mM dNTP (free Mg2+ = 0.7 mM)
    with_dntp = TMSettings(divalent_salt_conc=1.5, dntp_conc=0.8)
    # 50 mM Na+, 1.5 mM Mg2+, 2.0 mM dNTP (excess dNTP, free Mg2+ = 0.0 mM)
    excess_dntp = TMSettings(divalent_salt_conc=1.5, dntp_conc=2.0)
    # 50 mM Na+, 0.0 mM Mg2+, 0.0 mM dNTP (0 divalent)
    zero_mg = TMSettings(divalent_salt_conc=0.0, dntp_conc=0.0)

    seq = Primer("TAATACGACTCACTATAGGG")

    tm_no_dntp = calculate_tm_santalucia_1998_owczarzy_2008(seq, no_dntp)
    tm_with_dntp = calculate_tm_santalucia_1998_owczarzy_2008(seq, with_dntp)
    tm_excess = calculate_tm_santalucia_1998_owczarzy_2008(seq, excess_dntp)
    tm_zero_mg = calculate_tm_santalucia_1998_owczarzy_2008(seq, zero_mg)

    # dNTP chelation reduces available Mg2+, lowering Tm
    assert tm_with_dntp < tm_no_dntp
    # Excess dNTP eliminates free Mg2+, matching 0 divalent salt
    assert tm_excess == pytest.approx(tm_zero_mg)


def test_symmetry_correction() -> None:
    """Test that self-complementary sequences undergo symmetry correction."""
    # Self-complementary sequence: CGCGCGCGCGCG
    self_comp_seq = Primer("CGCGCGCGCGCG")
    tm = calculate_tm_santalucia_1998_owczarzy_2008(
        self_comp_seq, GLOBAL_TM_SETTINGS
    )
    assert tm == pytest.approx(66.0445, abs=0.5)


def test_negative_concentrations_safeguard() -> None:
    """Test negative salt/DNA concentrations do not crash calculations."""
    neg_settings = TMSettings(
        monovalent_salt_conc=-50.0,
        divalent_salt_conc=-1.5,
        dntp_conc=-0.2,
        dna_conc=-10.0,
    )
    tm_std = calculate_tm_santalucia_1998_owczarzy_2008(
        Primer("TAATACGACTCACTATAGGG"), neg_settings
    )
    assert tm_std > -273.15

    tm_amp4 = calculate_tm_lander_amplify4(primer_11bp, neg_settings)
    assert tm_amp4 > 0.0


def test_positive_mg_negative_dntp() -> None:
    """Test positive divalent salt with negative dNTP concentration."""
    pos_mg_neg_dntp = TMSettings(
        divalent_salt_conc=1.5,
        dntp_conc=-0.2,
    )
    pos_mg_zero_dntp = TMSettings(
        divalent_salt_conc=1.5,
        dntp_conc=0.0,
    )
    primer = Primer("TAATACGACTCACTATAGGG")
    tm_neg = calculate_tm_santalucia_1998_owczarzy_2008(primer, pos_mg_neg_dntp)
    tm_zero = calculate_tm_santalucia_1998_owczarzy_2008(
        primer, pos_mg_zero_dntp
    )
    assert tm_neg == pytest.approx(tm_zero)


def test_ambiguous_base_gc_content() -> None:
    """Test expected GC contribution for IUPAC degenerate bases S and W."""
    seq_s = Primer("ACGTACGTACS")
    seq_w = Primer("ACGTACGTACW")
    tm_s = calculate_tm_santalucia_1998_owczarzy_2008(seq_s, GLOBAL_TM_SETTINGS)
    tm_w = calculate_tm_santalucia_1998_owczarzy_2008(seq_w, GLOBAL_TM_SETTINGS)
    assert tm_s > tm_w


def test_calculate_tm_zero_enthalpy_and_denominators() -> None:
    """Test defensive error handling for zero thermodynamic denominators."""
    from unittest.mock import patch

    # 1. Zero dH in SantaLucia calculation (line 196)
    with patch.dict(
        "amplifyp.melting._NN_THERMO_DATA_TUPLE",
        {("A", "A"): (-920.0, -22.2)},
    ):
        with pytest.raises(
            InsufficientThermodynamicDataError,
            match="Invalid sequence: lacks standard thermodynamic base pairs",
        ):
            calculate_tm_santalucia_1998_owczarzy_2008(Primer("AAAAAA"))

    # 2. Zero calculated base Tm in Kelvin (line 202)
    import math

    orig_isclose = math.isclose

    call_count_tm = 0

    def isclose_base_tm(a: float, b: float, **kw: Any) -> bool:
        nonlocal call_count_tm
        call_count_tm += 1
        if call_count_tm == 4:
            return True
        return bool(orig_isclose(a, b, **kw))

    with patch("amplifyp.melting.math.isclose", side_effect=isclose_base_tm):
        with pytest.raises(
            InsufficientThermodynamicDataError,
            match="Calculated base Tm is zero",
        ):
            calculate_tm_santalucia_1998_owczarzy_2008(Primer("ACGTACGT"))

    # 3. Denominator with salt correction is zero (monovalent mode, line 215)
    s_mono = TMSettings(monovalent_salt_conc=50.0, divalent_salt_conc=0.0)
    call_count_mono = 0

    def isclose_mono(a: float, b: float, **kw: Any) -> bool:
        nonlocal call_count_mono
        call_count_mono += 1
        if call_count_mono == 6:
            return True
        return bool(orig_isclose(a, b, **kw))

    with patch("amplifyp.melting.math.isclose", side_effect=isclose_mono):
        with pytest.raises(
            InsufficientThermodynamicDataError,
            match="Denominator with salt correction is zero",
        ):
            calculate_tm_santalucia_1998_owczarzy_2008(
                Primer("ACGTACGT"), settings=s_mono
            )

    # 4. Denominator with salt correction is zero
    # (monovalent dominant mode, line 231)
    s_ratio = TMSettings(monovalent_salt_conc=100.0, divalent_salt_conc=0.0001)
    call_count_ratio = 0

    def isclose_ratio(a: float, b: float, **kw: Any) -> bool:
        nonlocal call_count_ratio
        call_count_ratio += 1
        if call_count_ratio == 6:
            return True
        return bool(orig_isclose(a, b, **kw))

    with patch("amplifyp.melting.math.isclose", side_effect=isclose_ratio):
        with pytest.raises(
            InsufficientThermodynamicDataError,
            match="Denominator with salt correction is zero",
        ):
            calculate_tm_santalucia_1998_owczarzy_2008(
                Primer("ACGTACGT"), settings=s_ratio
            )

    # 5. Inverse Tm with salt correction is zero (Mg mode, line 271)
    s_mg = TMSettings(monovalent_salt_conc=10.0, divalent_salt_conc=1.5)
    call_count_mg = 0

    def isclose_mg(a: float, b: float, **kw: Any) -> bool:
        nonlocal call_count_mg
        call_count_mg += 1
        if call_count_mg == 6:
            return True
        return bool(orig_isclose(a, b, **kw))

    with patch("amplifyp.melting.math.isclose", side_effect=isclose_mg):
        with pytest.raises(
            InsufficientThermodynamicDataError,
            match="Inverse Tm with salt correction is zero",
        ):
            calculate_tm_santalucia_1998_owczarzy_2008(
                Primer("ACGTACGT"), settings=s_mg
            )

    # 6. Amplify4 zero denominator (line 340)
    log_dna = 1.987 * math.log(50.0 / 4e9)

    target_e = 10.0 * log_dna - 108.0
    s_amp4 = TMSettings()
    s_amp4.entropy["C", "A"] = target_e
    with pytest.raises(
        InsufficientThermodynamicDataError,
        match="Denominator is zero in Tm calculation \\(Amplify4\\)",
    ):
        calculate_tm_lander_amplify4(Primer("AC"), settings=s_amp4)


def test_calculate_tm_dispatcher_and_fallback() -> None:
    """Test calculate_tm dispatch across methods and Amplify4 fallback."""
    std_primer = primer_1701
    assert calculate_tm(
        std_primer, method=TmMethod.SANTALUCIA_1998_OWCZARZY_2008
    ) == pytest.approx(
        calculate_tm_santalucia_1998_owczarzy_2008(
            std_primer, GLOBAL_TM_SETTINGS
        )
    )
    assert calculate_tm(
        std_primer, method=TmMethod.LANDER_AMPLIFY4
    ) == pytest.approx(
        calculate_tm_lander_amplify4(std_primer, GLOBAL_TM_SETTINGS)
    )

    deg_primer = Primer("ACGTNNACGT")
    expected_sl = calculate_tm_santalucia_1998_owczarzy_2008(
        deg_primer, GLOBAL_TM_SETTINGS
    )
    assert calculate_tm(deg_primer, fallback_to_amplify4=True) == pytest.approx(
        expected_sl
    )
    assert calculate_tm(
        deg_primer, fallback_to_amplify4=False
    ) == pytest.approx(expected_sl)
    assert calculate_tm(
        deg_primer, method=TmMethod.LANDER_AMPLIFY4
    ) == pytest.approx(
        calculate_tm_lander_amplify4(deg_primer, GLOBAL_TM_SETTINGS)
    )

    with pytest.raises(ValueError):
        calculate_tm(std_primer, method="Unknown Method")

    # Fallback path when SantaLucia/Owczarzy raises for unsupported bases
    with patch(
        "amplifyp.melting.calculate_tm_santalucia_1998_owczarzy_2008",
        side_effect=InsufficientThermodynamicDataError("unsupported"),
    ):
        assert calculate_tm(
            deg_primer, fallback_to_amplify4=True
        ) == pytest.approx(
            calculate_tm_lander_amplify4(deg_primer, GLOBAL_TM_SETTINGS)
        )
        with pytest.raises(InsufficientThermodynamicDataError):
            calculate_tm(deg_primer, fallback_to_amplify4=False)

    # Degenerate IUPAC bases other than N are supported by the default model
    # while still raising InsufficientThermodynamicDataError in Amplify4.
    for deg_seq in ("ACGTRACGT", "ACGTYACGT", "ACGTMACGT", "ACGTVACGT"):
        iupac_primer = Primer(deg_seq)
        assert calculate_tm(iupac_primer) == pytest.approx(
            calculate_tm_santalucia_1998_owczarzy_2008(
                iupac_primer, GLOBAL_TM_SETTINGS
            )
        )
        with pytest.raises(InsufficientThermodynamicDataError):
            calculate_tm_lander_amplify4(iupac_primer, GLOBAL_TM_SETTINGS)
        with pytest.raises(InsufficientThermodynamicDataError):
            calculate_tm(iupac_primer, method=TmMethod.LANDER_AMPLIFY4)


def test_all_iupac_bases_and_equimolar_concentration() -> None:
    """Test all 15 IUPAC codes and equimolar concentration division."""
    all_iupac = Primer("".join(IUPAC_BASES.keys()))
    tm_all = calculate_tm_santalucia_1998_owczarzy_2008(
        all_iupac, GLOBAL_TM_SETTINGS
    )
    assert tm_all > -273.15

    # Equimolar concentration check: ACGTNACGT (fold=4) has variants at
    # dna_conc / 4, and S/W ordering holds for GC- vs AT-rich variants.
    p_n = Primer("ACGTNACGTACGT")
    per_species_settings = replace(
        GLOBAL_TM_SETTINGS, dna_conc=GLOBAL_TM_SETTINGS.dna_conc / 4.0
    )
    variant_tms = [
        calculate_tm_santalucia_1998_owczarzy_2008(
            Primer(f"ACGT{base}ACGTACGT"), per_species_settings
        )
        for base in ("A", "C", "G", "T")
    ]
    tm_range = calculate_tm_range(p_n, GLOBAL_TM_SETTINGS)
    assert tm_range.variant_count == 4
    assert tm_range.is_approximate is False
    assert tm_range.min_tm == pytest.approx(min(variant_tms))
    assert tm_range.max_tm == pytest.approx(max(variant_tms))
    assert tm_range.mean_tm == pytest.approx(
        sum(variant_tms) / len(variant_tms)
    )
    mean_approx = calculate_tm_santalucia_1998_owczarzy_2008(
        p_n, GLOBAL_TM_SETTINGS
    )
    assert abs(mean_approx - tm_range.mean_tm) < 0.5


def test_calculate_tm_range_enumeration_cap_and_dp_bounds() -> None:
    """Test exact enumeration up to cap and O(L) DP bounds above cap."""
    # Exact at boundary (4**3 = 64 == max_variants)
    p64 = Primer("ACGTNNNACGTA")
    exact_at_cap = calculate_tm_range(p64, GLOBAL_TM_SETTINGS, max_variants=64)
    assert exact_at_cap.variant_count == 64
    assert exact_at_cap.is_approximate is False

    # Approximate above boundary (64 > 16)
    approx = calculate_tm_range(p64, GLOBAL_TM_SETTINGS, max_variants=16)
    assert approx.variant_count == 64
    assert approx.is_approximate is True
    assert approx.min_tm <= approx.mean_tm <= approx.max_tm
    assert approx.min_tm == pytest.approx(exact_at_cap.min_tm, abs=0.5)
    assert approx.max_tm == pytest.approx(exact_at_cap.max_tm, abs=0.5)

    # Very high redundancy avoids float overflow (4**600 > 1e361)
    ultra = calculate_tm_range(Primer("N" * 600), GLOBAL_TM_SETTINGS)
    assert isinstance(ultra, TmRange)
    assert ultra.is_approximate is True
    assert ultra.min_tm < ultra.mean_tm < ultra.max_tm


def test_calculate_tm_range_edge_cases_and_methods() -> None:
    """Test calculate_tm_range on short, canonical, Lander, and bad inputs."""
    assert calculate_tm_range(Primer("")) == TmRange(
        min_tm=0.0,
        mean_tm=0.0,
        max_tm=0.0,
        is_approximate=False,
        variant_count=1,
    )
    assert calculate_tm_range(Primer("N")) == TmRange(
        min_tm=0.0,
        mean_tm=0.0,
        max_tm=0.0,
        is_approximate=False,
        variant_count=4,
    )

    std_tm = calculate_tm_santalucia_1998_owczarzy_2008(
        primer_1701, GLOBAL_TM_SETTINGS
    )
    assert calculate_tm_range(primer_1701, GLOBAL_TM_SETTINGS) == TmRange(
        min_tm=std_tm,
        mean_tm=std_tm,
        max_tm=std_tm,
        is_approximate=False,
        variant_count=1,
    )

    amp4_tm = calculate_tm_lander_amplify4(
        Primer("ACGTNNACGT"), GLOBAL_TM_SETTINGS
    )
    assert calculate_tm_range(
        Primer("ACGTNNACGT"),
        GLOBAL_TM_SETTINGS,
        method=TmMethod.LANDER_AMPLIFY4,
    ) == TmRange(
        min_tm=amp4_tm,
        mean_tm=amp4_tm,
        max_tm=amp4_tm,
        is_approximate=False,
        variant_count=1,
    )

    with pytest.raises(ValueError, match="max_variants must be at least 1"):
        calculate_tm_range(primer_1701, max_variants=0)

    bad_primer = MagicMock(spec=Primer)
    bad_primer.seq_upper = "AC-T"
    with pytest.raises(InsufficientThermodynamicDataError):
        calculate_tm_range(bad_primer, fallback_to_amplify4=False)
    with patch(
        "amplifyp.melting.calculate_tm_lander_amplify4", return_value=42.0
    ):
        assert calculate_tm_range(
            bad_primer, fallback_to_amplify4=True
        ) == TmRange(
            min_tm=42.0,
            mean_tm=42.0,
            max_tm=42.0,
            is_approximate=False,
            variant_count=1,
        )
