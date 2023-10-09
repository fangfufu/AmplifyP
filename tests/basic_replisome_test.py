# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNAType
from amplifyp.replisome import Replisome


def test_replisome_simple() -> None:
    """A simple replisome test."""
    target = DNA("GTGTACTCAGTTCCATAAACGAGCTATTAGATATGAGGTCCGTAGATTGAAAAGGGTGA")
    primer = DNA("GTGTACT", DNAType.PRIMER)
    min_overlap = 3
    replisome = Replisome(target, primer, min_overlap=min_overlap)
    test_index = replisome.target_index_limit.stop - len(primer) + min_overlap
    test_str = replisome.target[replisome.replicon_slice(test_index)]
    assert test_str == "AGTGGGA"
