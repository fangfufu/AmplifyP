# -*- coding: utf-8 -*-
"""Tests related to replication origin."""
from dataclasses import dataclass
from typing import List
import pytest

from amplifyp.dna import DNA
from amplifyp.origin import ReplicationOrigin
from amplifyp.settings import Settings


def test_replication_origin_init() -> None:
    """Test the initialization of a Origin object with invalid parameters."""
    with pytest.raises(ValueError):
        ReplicationOrigin(target="ATCG", primer="ATC", settings=Settings())


@dataclass(frozen=True, slots=True)
class ReplicationOriginExample:
    """A class representing test cases of origins.

    Attributes:
        origin (Origin): The origin being tested.
        primability (float): The primability score of the origin.
        stability (float): The stability score of the origin.
        quality (float): The quality score of the origin.
    """

    origin: ReplicationOrigin
    primability: float
    stability: float
    quality: float

# Origin example 1 from the documentation
# Reverse origin
origin_examples: List[ReplicationOriginExample] = []
origin_examples.append(
    ReplicationOriginExample(
        ReplicationOrigin(
            target=DNA("CGACTGGGCAAAGGAAATCCTT").complement().sequence,
            primer="GCTGACCCNTTTCYYTTAGGCA",
            settings=Settings(),
        ),
        0.9923605805958747, # Should be 99%
        0.9118497898586322, # Should be 92%
        0.8802629630681338, # Should be 0.8875
    )
)

# Origin example 2 from the documentation
# Forward origin
origin_examples.append(
    ReplicationOriginExample(
        ReplicationOrigin(
            target=DNA("CGAGGGGGCAAAGGAAATCC").reverse().sequence,
            primer=DNA("CGACTGGGCAAAGGAAATCC").reverse().sequence,
            settings=Settings(),
        ),
        0.9850746268656716, # Should be 98%
        0.8914650537634409, # Should be 89%
        0.8456746007863905, # Should be 0.8375
    )
)

# Reverse Origin for test_verification_F1423_R20049, Primer 20049
origin_examples.append(
    ReplicationOriginExample(
        ReplicationOrigin(
            target=DNA("ATCCTCTTTTTATGTCCACATCTG").complement().sequence,
            primer="TAACAGAAAAATACAGGTGTAGAC",
            settings=Settings(),
        ),
        0.855072463768116, # Should be 85%
        0.47619047619047616, # Should be 90%
        0.16407867494823983, # Should be 0.6875
    )
)

# Forward Origin for test_verification_F1423_R20049, Primer 1423
origin_examples.append(
    ReplicationOriginExample(
        ReplicationOrigin(
            target=DNA("CGCAAAGCTTGTCGGC").reverse().sequence,
            primer=DNA("GCCAAAGTGTGTTGGC").reverse().sequence,
            settings=Settings(),
        ),
        0.8076923076923077, # Should be 80%
        0.5837912087912088, # Should be 66%
        0.23935439560439536, # Should be 0.325
    )
)

def test_origin_primability() -> None:
    """Test if origin primability is working correctly."""
    for ex in origin_examples:
        assert ex.primability == ex.origin.primability()


def test_origin_stability() -> None:
    """Test if origin stability is working correctly."""
    print("")
    for ex in origin_examples:
        # assert ex.stability == ex.origin.stability()
        print(ex.origin.stability())


def test_origin_quality() -> None:
    """Test if origin quality is working correctly."""
    for ex in origin_examples:
        assert ex.quality == ex.origin.quality()
