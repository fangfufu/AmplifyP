# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNAType
from amplifyp.replisome import Replisome


def test_replisome_replicon_slice_limit_stop() -> None:
    """Test replicon slicing at the end of the index limit."""
    target: DNA = DNA("GTGTACTCAGTTCCATAAACGAGCTATTAGATATGAGGTCCGTAGATTGAAAAGGGTGA")
    primer: DNA = DNA("GTGTACT", DNAType.PRIMER)
    target_str: str = target[0 : len(primer)].sequence[::-1] + len(primer) * "-"
    assert isinstance(target_str, str)
    for i in range(0, len(primer)):
        replisome = Replisome(target, primer, min_overlap=i)
        test_index = replisome.target_index_limit.stop - i
        test_str: str = replisome.target[replisome.replicon_slice(test_index)].sequence
        assert isinstance(test_str, str)
        assert test_str == target_str[len(primer) - i : 2 * len(primer) - i]
