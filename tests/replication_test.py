# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNAType, Primer
from amplifyp.replication import Repliconf
from amplifyp.settings import Settings


def test_repliconf_idx_to_template_idx_linear() -> None:
    """Test case for the `idx_repliconf_to_template` for linear template."""
    # Create an instance of the class
    target = DNA("ATCGATCG")
    primer = Primer("CGAT")
    min_overlap = 2
    config = Settings(min_overlap=min_overlap)

    replication = Repliconf(target, primer, config)

    origin_idx = 3
    expected_result = 4
    assert replication.idx_repliconf_to_template(origin_idx) == expected_result


def test_repliconf_idx_to_template_idx_circular() -> None:
    """Test case for the `idx_repliconf_to_template` for circular template."""
    target = DNA("ATCGATCG", DNAType.CIRCULAR)
    primer = Primer("CGAT")
    min_overlap = 2
    config = Settings(min_overlap=min_overlap)

    replication = Repliconf(target, primer, config)

    origin_idx = 7
    expected_result = 2
    assert replication.idx_repliconf_to_template(origin_idx) == expected_result
