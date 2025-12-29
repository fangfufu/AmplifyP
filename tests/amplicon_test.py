# -*- coding: utf-8 -*-
"""Tests amplicon generation."""

import pytest
from amplifyp.amplicon import AmpliconGenerator
from amplifyp.dna import DNA, DNAType, Primer
from amplifyp.repliconf import Repliconf
from amplifyp.settings import DEFAULT_SETTINGS


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
