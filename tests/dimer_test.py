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

import pytest

from amplifyp.dimer import PrimerDimer, PrimerDimerGenerator
from amplifyp.dna import Primer
from amplifyp.settings import PrimerDimerSettings
from tests.examples.amplify4_examples import (
    primer_11bp,
    primer_1701,
    primer_10289,
    primer_10290,
)


def test_primer_dimer_generator() -> None:
    """Test PrimerDimerGenerator using real examples from Amplify4."""
    primer_dimer_generator = PrimerDimerGenerator()

    # Real examples from Amplify4
    result_10289_10289 = primer_dimer_generator.generate_primer_dimer(
        primer_10289, primer_10289
    )
    assert result_10289_10289.quality == 160
    assert result_10289_10289.overlap == 10

    result_10290_10290 = primer_dimer_generator.generate_primer_dimer(
        primer_10290, primer_10290
    )
    assert result_10290_10290.quality == 200
    assert result_10290_10290.overlap == 12

    result_10290_11bp = primer_dimer_generator.generate_primer_dimer(
        primer_10290, primer_11bp
    )
    assert result_10290_11bp.quality == 70
    assert result_10290_11bp.overlap == 6

    result_1701_1701 = primer_dimer_generator.generate_primer_dimer(
        primer_1701, primer_1701
    )
    assert result_1701_1701.quality == 120
    assert result_1701_1701.overlap == 8

    result_11bp_11bp = primer_dimer_generator.generate_primer_dimer(
        primer_11bp, primer_11bp
    )
    assert result_11bp_11bp.quality == 60
    assert result_11bp_11bp.overlap == 9

    primer_dimer_generator.add_primer(primer_10289)
    primer_dimer_generator.add_primer(primer_10290)
    primer_dimer_generator.add_primer(primer_1701)
    primer_dimer_generator.add_primer(primer_11bp)
    primer_dimer_generator.analyse_primers()
    assert len(primer_dimer_generator.primer_dimers) == 4
    assert result_10289_10289 in primer_dimer_generator.primer_dimers
    assert result_10290_10290 in primer_dimer_generator.primer_dimers
    assert result_10290_11bp in primer_dimer_generator.primer_dimers
    assert result_1701_1701 in primer_dimer_generator.primer_dimers
    assert result_11bp_11bp not in primer_dimer_generator.primer_dimers


def test_generator_management() -> None:
    """Test adding, removing, clearing primers and analysis state."""
    generator = PrimerDimerGenerator()
    p1 = primer_10289
    p2 = primer_10290

    # Test add
    generator.add_primer(p1)
    generator.add_primer(p2)
    assert len(generator.primers) == 2
    assert p1 in generator.primers
    assert p2 in generator.primers
    assert not generator.analysed

    # Test remove
    generator.remove_primer(p1)
    assert len(generator.primers) == 1
    assert p1 not in generator.primers
    assert p2 in generator.primers

    # Test clear
    generator.add_primer(p1)
    generator.analyse_primers()
    assert generator.analysed
    assert (
        len(generator.primer_dimers) > 0
    )  # Should have self-dimers if they pass threshold

    generator.add_primer(p2)
    assert not generator.analysed

    generator.analyse_primers()
    assert generator.analysed
    generator.remove_primer(p2)
    assert not generator.analysed

    generator.clear()
    assert len(generator.primers) == 0
    assert len(generator.primer_dimers) == 0
    assert not generator.analysed


def test_equal_length_primer_symmetry() -> None:
    """Test that generate_primer_dimer is symmetric for equal length primers."""
    p_a = Primer("TTTTTCCCCC", name="p_a")
    p_b = Primer("GGGGGAAAAA", name="p_b")
    generator = PrimerDimerGenerator()
    dimer1 = generator.generate_primer_dimer(p_a, p_b)
    dimer2 = generator.generate_primer_dimer(p_b, p_a)

    assert dimer1.quality == dimer2.quality
    assert dimer1.overlap == dimer2.overlap


def test_custom_settings() -> None:
    """Test using custom settings for the generator."""
    # Create settings with a very high threshold so nothing should pass
    high_threshold_settings = PrimerDimerSettings(threshold=1000.0)
    generator = PrimerDimerGenerator(settings=high_threshold_settings)

    generator.add_primer(primer_10289)
    generator.add_primer(primer_10290)
    generator.analyse_primers()

    assert len(generator.primer_dimers) == 0

    # Create settings with low overlap requirement
    # primer_11bp self dimer has Q=60, overlap=9. Default threshold=60.
    # Let's set threshold higher than 60 to filter it out.
    # Default settings: threshold 60. result_11bp_11bp Q=60.
    # Logic is `> threshold` (strictly greater). So 60 > 60 is False.
    # That explains why it wasn't in the original test:
    # `assert result_11bp_11bp not in primer_dimer_generator.primer_dimers`

    # Let's lower threshold to 59
    low_threshold_settings = PrimerDimerSettings(threshold=59.0)
    generator_low = PrimerDimerGenerator(settings=low_threshold_settings)
    generator_low.add_primer(primer_11bp)
    generator_low.analyse_primers()

    # Should now be included
    # We need to construct what the result would be or check existence
    assert len(generator_low.primer_dimers) == 1
    assert generator_low.primer_dimers[0].primer_1 == primer_11bp
    assert generator_low.primer_dimers[0].primer_2 == primer_11bp


def test_edge_cases() -> None:
    """Test edge cases like short primers, no overlap potential, etc."""
    generator = PrimerDimerGenerator()

    # Very short primer - unlikely to have high score
    short_p = Primer("ATCG", name="short")

    # Compare with itself
    res = generator.generate_primer_dimer(short_p, short_p)
    # 4 bases. Max score roughly 4*30 = 120? But alignments might be poor.
    # Just check it doesn't crash
    assert isinstance(res, PrimerDimer)
    assert res.primer_1 == short_p
    assert res.primer_2 == short_p

    # Primers with no complementarity (poly-A vs poly-A)
    # Default weights: A-A mismatches are penalty (-20 or similar)
    pA = Primer("AAAAA", name="polyA")
    res_A = generator.generate_primer_dimer(pA, pA)
    # Should be low quality
    assert res_A.quality < 0


def test_sorting_and_filtering() -> None:
    """Test that analyse_primers filters and sorts correctly."""
    # We'll use 3 primers that produce dimers of varying quality
    # p1-p1: Q=100 (hypothetical)
    # p2-p2: Q=80
    # p3-p3: Q=50
    # Threshold=60.

    # We can use the real examples since we know their scores.
    # 10290-10290: 200
    # 10289-10289: 160
    # 1701-1701: 120
    # 11bp-11bp: 60 (filtered out by >60)

    generator = PrimerDimerGenerator()
    generator.add_primer(primer_10290)
    generator.add_primer(primer_10289)
    generator.add_primer(primer_1701)
    generator.add_primer(primer_11bp)

    generator.analyse_primers()

    results = generator.primer_dimers

    # Verify sorting (descending quality)
    qualities = [pd.quality for pd in results]
    assert qualities == sorted(qualities, reverse=True)

    # Verify filtering
    # 11bp self dimer has Q=60. Threshold is 60. > check excludes it.
    for pd in results:
        assert pd.quality > 60.0
        assert pd.overlap > 3  # default min overlap


def test_primer_dimer_attributes() -> None:
    """Test PrimerDimer dataclass attributes."""
    pd = PrimerDimer(
        primer_1=primer_11bp,
        primer_2=primer_1701,
        quality=123.4,
        overlap=5,
        p1_pos=2,
    )
    assert pd.primer_1 == primer_11bp
    assert pd.primer_2 == primer_1701
    assert pd.quality == pytest.approx(123.4)
    assert pd.overlap == 5
    assert pd.p1_pos == 2


def test_primer_order_swap() -> None:
    """Test reorder parameter in generate_primer_dimer."""
    generator = PrimerDimerGenerator()
    p_long = Primer("AAAAATTTTT", name="long")
    p_short = Primer("AAAAA", name="short")

    # Default reorder=True reorders so shorter primer is primer_1
    res_default = generator.generate_primer_dimer(p_long, p_short)
    assert res_default.primer_1 == p_short
    assert res_default.primer_2 == p_long

    # reorder=False preserves primer order
    res_no_reorder = generator.generate_primer_dimer(
        p_long, p_short, reorder=False
    )
    assert res_no_reorder.primer_1 == p_long
    assert res_no_reorder.primer_2 == p_short
    assert res_default.quality == res_no_reorder.quality
    assert res_default.overlap == res_no_reorder.overlap
    assert res_no_reorder.binding_strength_str == "|||||"


@pytest.mark.parametrize(  # type: ignore[untyped-decorator]
    (
        "p1, p2, expected_overlap, expected_quality, "
        "expected_p1_pos, expected_strength"
    ),
    [
        (primer_10290, primer_10290, 12, 200.0, 8, "| |||||||| |"),
        (primer_10290, primer_11bp, 6, 70.0, 14, "|  |||"),
        (primer_10290, primer_10289, 7, -10.0, 13, " | |  |"),
        (primer_10290, primer_1701, 2, 60.0, 18, "||"),
        (primer_11bp, primer_11bp, 9, 60.0, 2, "|||   |||"),
        (primer_11bp, primer_10289, 11, 30.0, 7, " ||  ||  ||"),
        (primer_11bp, primer_1701, 11, 10.0, 2, "   |||  || "),
        (primer_10289, primer_10289, 10, 160.0, 10, "|| |||| ||"),
        (primer_10289, primer_1701, 6, 50.0, 14, "| ||| "),
        (primer_1701, primer_1701, 8, 120.0, 12, "|||  |||"),
    ],
)
def test_primer_dimer_binding_strength(
    p1: Primer,
    p2: Primer,
    expected_overlap: int,
    expected_quality: float,
    expected_p1_pos: int,
    expected_strength: str,
) -> None:
    """Test binding strength string generation and dimer calculations for
    specified primers.
    """
    generator = PrimerDimerGenerator()
    dimer = generator.generate_primer_dimer(p1, p2)

    assert dimer.overlap == expected_overlap
    assert dimer.quality == expected_quality
    assert dimer.p1_pos == expected_p1_pos
    assert dimer.binding_strength_str == expected_strength


def test_primer_dimer_bonds_key_error_and_colon() -> None:
    """Test bond symbol generation with unknown bases and custom threshold."""
    from amplifyp.settings import BasePairWeightsTbl

    # 4x4 table without degenerate base 'M' causes KeyError
    # in weights -> score 0.0 -> ':' symbol
    tbl = BasePairWeightsTbl(
        row="ACGT",
        col="ACGT",
        weight=[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
    )
    s_missing = PrimerDimerSettings(weights=tbl)
    pd_unknown = PrimerDimer(
        primer_1=Primer("MMM"),
        primer_2=Primer("TTT"),
        quality=10.0,
        overlap=3,
        p1_pos=0,
        settings=s_missing,
    )
    assert pd_unknown.binding_strength_str == ":::"

    # Custom settings with negative score results in space ' '
    tbl_neg = BasePairWeightsTbl(
        row="ACGT",
        col="ACGT",
        weight=[
            [-10, -10, -10, -10],
            [-10, -10, -10, -10],
            [-10, -10, -10, -10],
            [-10, -10, -10, -10],
        ],
    )
    s_neg = PrimerDimerSettings(weights=tbl_neg)
    pd_space = PrimerDimer(
        primer_1=Primer("AAA"),
        primer_2=Primer("TTT"),
        quality=-10.0,
        overlap=3,
        p1_pos=0,
        settings=s_neg,
    )
    assert pd_space.binding_strength_str == "   "


def test_calculate_dimer_stats_key_error() -> None:
    """Test calculate_dimer_stats with non-standard nucleotide characters."""
    from amplifyp.dimer import calculate_dimer_stats

    # 'Z' is not in standard primer dimer weights
    best_q, best_pos, overlap = calculate_dimer_stats(
        "AAAZAAA", "TTTTTTT", 7, 7
    )
    assert isinstance(best_q, float)
    assert isinstance(best_pos, int)
    assert isinstance(overlap, int)


def test_primer_dimer_reverse_binding_swap() -> None:
    """Test symmetric primers where reverse orientation has higher score."""
    custom_settings = PrimerDimerSettings()

    custom_settings.weights["A", "T"] = 10.0
    custom_settings.weights["T", "A"] = 25.0

    p1 = Primer("AAAAA", name="P1")
    p2 = Primer("TTTTT", name="P2")

    generator = PrimerDimerGenerator(settings=custom_settings)
    res = generator.generate_primer_dimer(p1, p2)
    assert res.primer_1 == p2
    assert res.primer_2 == p1

    generator.add_primer(p1)
    generator.add_primer(p2)
    generator.analyse_primers()
    assert len(generator.primer_dimers) >= 1
    found_dimer = generator.primer_dimers[0]
    assert found_dimer.primer_1 == p2
    assert found_dimer.primer_2 == p1
