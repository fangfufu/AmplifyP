"""Tests amplicon generation."""

import pytest

from amplifyp.amplicon import AmpliconGenerator
from amplifyp.dna import DNA, DNAType, Primer
from amplifyp.repliconf import Repliconf
from amplifyp.settings import DEFAULT_SETTINGS, Settings


def test_generator_init() -> None:
    """Test initialization of AmpliconGenerator."""
    template = DNA("ATGC", DNAType.LINEAR, "Template 1")
    generator = AmpliconGenerator(template)
    assert generator.template == template
    assert not generator.repliconfs


def test_generator_add_success() -> None:
    """Test adding a compatible Repliconf."""
    template = DNA("ATGC", DNAType.LINEAR, "Template 1")
    primer = Primer("AT", "Primer 1")
    repliconf = Repliconf(template, primer, DEFAULT_SETTINGS)

    generator = AmpliconGenerator(template)
    generator.add(repliconf)

    assert len(generator.repliconfs) == 1
    assert generator.repliconfs[0] == repliconf


def test_generator_add_failure() -> None:
    """Test adding a Repliconf with a different template raises ValueError."""
    template1 = DNA("ATGC", DNAType.LINEAR, "Template 1")
    template2 = DNA("CGTA", DNAType.LINEAR, "Template 2")
    primer = Primer("AT", "Primer 1")
    repliconf = Repliconf(template2, primer, DEFAULT_SETTINGS)

    generator = AmpliconGenerator(template1)

    with pytest.raises(ValueError, match="The Repliconf contains a different template"):
        generator.add(repliconf)


def test_generate_amplicons() -> None:
    """Test generating amplicons with multiple configurations."""
    # Template: AAAAA + (GT)*10 + GGGGG (30 bp)
    seq = "AAAAA" + "GT" * 10 + "GGGGG"
    template = DNA(seq, DNAType.LINEAR, "CleanTemplate")

    primer_fwd = Primer("AAAAA", "PF")
    primer_rev = Primer("CCCCC", "PR")

    settings = Settings()
    settings.primability_cutoff = 0.0
    settings.stability_cutoff = -100.0

    rc_fwd = Repliconf(template, primer_fwd, settings)
    rc_fwd.search()

    rc_rev = Repliconf(template, primer_rev, settings)
    rc_rev.search()

    generator = AmpliconGenerator(template)
    generator.add(rc_fwd)
    generator.add(rc_rev)

    amplicons = generator.generate_amplicons()
    assert len(amplicons) >= 1

    longest_amp = max(amplicons, key=lambda a: len(a.sequence))
    assert len(longest_amp.sequence) == 30
    assert longest_amp.sequence.seq == seq
