# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from dataclasses import dataclass
from typing import List
import pytest

from amplifyp.dna import DNA, DNAType
from amplifyp.replisome import Replisome, Replicon, ReplicationConfig


def test_replisome_linear_replicon_slice_limit_stop() -> None:
    """Test replicon slicing for linear DNA at the end of the index limit."""
    target: DNA = DNA("GTGTACTCAGTTCCATAAACGAGCTATTAGATATGAGGTCCGTAGATTGAAAAGGGTGA")
    primer: DNA = DNA("GTGTACT", DNAType.PRIMER)
    target_str: str = target[0 : len(primer)].sequence[::-1] + len(primer) * "-"
    assert isinstance(target_str, str)
    for i in range(0, len(primer)):
        replication_config = ReplicationConfig(min_overlap=i)
        replisome = Replisome(target, primer, replication_config)
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
        replication_config = ReplicationConfig(min_overlap=i)
        replisome = Replisome(target, primer, replication_config)
        test_index = replisome.range_limit.start
        test_str: str = replisome.target[replisome.replicon_slice(test_index)].sequence
        assert isinstance(test_str, str)
        assert test_str == target_str[i : len(primer) + i]


def test_replicon_init() -> None:
    """Test the initialization of a Replicon object with invalid parameters."""
    with pytest.raises(ValueError):
        Replicon(target="ATCG", primer="ATC", replication_config=ReplicationConfig())


@dataclass
class SingleRepliconResults:
    """A class representing test cases of replicons.

    Attributes:
        replicon (Replicon): The replicon being tested.
        primability (float): The primability score of the replicon.
        stability (float): The stability score of the replicon.
        quality (float): The quality score of the replicon.
    """

    replicon: Replicon
    primability: float
    stability: float
    quality: float


single_replicon_results: List[SingleRepliconResults] = []
single_replicon_results.append(
    SingleRepliconResults(
        Replicon(
            target=DNA("CGACTGGGCAAAGGAAATCCTT").complement().reverse().sequence,
            primer="GCTGACCCNTTTCYYTTAGGCA",
            replication_config=ReplicationConfig(),
        ),
        0.9923605805958747,
        0.9118497898586322,
        0.8802629630681338,
    )
)

single_replicon_results.append(
    SingleRepliconResults(
        Replicon(
            target=DNA("CGAGGGGGCAAAGGAAATCC").reverse().sequence,
            primer=DNA("CGACTGGGCAAAGGAAATCC").reverse().sequence,
            replication_config=ReplicationConfig(),
        ),
        0.9850746268656716,
        0.8914650537634409,
        0.8456746007863905,
    )
)


def test_replicon_primability() -> None:
    """Test if replicon primability is working correctly."""
    for i in single_replicon_results:
        assert i.primability == i.replicon.primability


def test_replicon_stability() -> None:
    """Test if replicon stability is working correctly."""
    for i in single_replicon_results:
        assert i.stability == i.replicon.stability


def test_replicon_quality() -> None:
    """Test if replicon quality is working correctly."""
    for i in single_replicon_results:
        assert i.quality == i.replicon.quality
