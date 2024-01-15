# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNAType, Primer
from amplifyp.replication import Repliconf
from amplifyp.settings import Settings


def test_repliconf_idx_to_template_idx_linear() -> None:
    """Test case for the `idx_repliconf_to_template` for linear template."""
    # Create an instance of the class
    template = DNA("ATCGATCG")
    primer = Primer("CGAT")
    # Target: 3'ATCGATCG----5'
    # Primer:   ##CGAT######
    config = Settings()

    replication = Repliconf(template, primer, config)

    repliconf_idx = 6
    template_idx = 2
    assert replication.idx_repliconf_to_template(repliconf_idx) == template_idx


def test_repliconf_idx_to_template_idx_circular() -> None:
    """Test case for the `idx_repliconf_to_template` for circular template."""
    template = DNA("ATCGATCG", DNAType.CIRCULAR)
    primer = Primer("CGAT")
    # Target: 3'ATCGATCGATCG5'
    # Primer:   AT####CG####
    config = Settings()

    replication = Repliconf(template, primer, config)

    repliconf_idx = 10
    template_idx = -2
    assert replication.idx_repliconf_to_template(repliconf_idx) == template_idx
