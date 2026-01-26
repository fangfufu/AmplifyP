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

    primer_dimer_generator.add(primer_10289)
    primer_dimer_generator.add(primer_10290)
    primer_dimer_generator.add(primer_1701)
    primer_dimer_generator.add(primer_11bp)
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
    generator.add(p1)
    generator.add(p2)
    assert len(generator.primers) == 2
    assert p1 in generator.primers
    assert p2 in generator.primers
    assert not generator.analysed

    # Test remove
    generator.remove(p1)
    assert len(generator.primers) == 1
    assert p1 not in generator.primers
    assert p2 in generator.primers

    # Test clear
    generator.add(p1)
    generator.analyse_primers()
    assert generator.analysed
    assert (
        len(generator.primer_dimers) > 0
    )  # Should have self-dimers if they pass threshold

    generator.clear()
    assert len(generator.primers) == 0
    assert len(generator.primer_dimers) == 0
    assert not generator.analysed


def test_custom_settings() -> None:
    """Test using custom settings for the generator."""
    # Create settings with a very high threshold so nothing should pass
    high_threshold_settings = PrimerDimerSettings(threshold=1000.0)
    generator = PrimerDimerGenerator(settings=high_threshold_settings)

    generator.add(primer_10289)
    generator.add(primer_10290)
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
    generator_low.add(primer_11bp)
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
    short_p = Primer("ATCG", "short")

    # Compare with itself
    res = generator.generate_primer_dimer(short_p, short_p)
    # 4 bases. Max score roughly 4*30 = 120? But alignments might be poor.
    # Just check it doesn't crash
    assert isinstance(res, PrimerDimer)
    assert res.primer_1 == short_p
    assert res.primer_2 == short_p

    # Primers with no complementarity (poly-A vs poly-A)
    # Default weights: A-A mismatches are penalty (-20 or similar)
    pA = Primer("AAAAA", "polyA")
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
    generator.add(primer_10290)
    generator.add(primer_10289)
    generator.add(primer_1701)
    generator.add(primer_11bp)

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
