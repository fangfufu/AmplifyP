# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNAType
from amplifyp.replisome import Replisome, ReplisomeConfig


def test_replisome_linear_replicon_slice_limit_stop() -> None:
    """Test replicon slicing for linear DNA at the end of the index limit."""
    target: DNA = DNA("GTGTACTCAGTTCCATAAACGAGCTATTAGATATGAGGTCCGTAGATTGAAAAGGGTGA")
    primer: DNA = DNA("GTGTACT", DNAType.PRIMER)
    target_str: str = target[0 : len(primer)].sequence[::-1] + len(primer) * "-"
    assert isinstance(target_str, str)
    for i in range(0, len(primer)):
        replisome_config = ReplisomeConfig(min_overlap=i)
        replisome = Replisome(target, primer, replisome_config)
        test_index = replisome.range_limit.stop - i
        test_str: str = replisome.target[replisome.replicon_slice(test_index)].sequence
        assert isinstance(test_str, str)
        assert test_str == target_str[len(primer) - i : 2 * len(primer) - i]


def test_replisome_linear_replicon_slice_limit_start() -> None:
    """Test replicon slicing for linear DNA at the start of the index limit."""
    target: DNA = DNA("ATTGTGCGATCCCTGCACCTCAGCTAAGGTAGCTACCAATATTTAGTTTCTAAGCCTTGC")
    primer: DNA = DNA("ATTGTGCGAT", DNAType.PRIMER)
    target_str: str = len(primer) * "-" + target.sequence[-len(primer) : :][::-1]
    assert isinstance(target_str, str)
    for i in range(0, len(primer)):
        replisome_config = ReplisomeConfig(min_overlap=i)
        replisome = Replisome(target, primer, replisome_config)
        test_index = replisome.range_limit.start
        test_str: str = replisome.target[replisome.replicon_slice(test_index)].sequence
        assert isinstance(test_str, str)
        assert test_str == target_str[i : len(primer) + i]
