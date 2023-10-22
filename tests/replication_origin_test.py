# -*- coding: utf-8 -*-
"""Tests related to replication origin."""
from dataclasses import dataclass
from typing import List
import pytest

from amplifyp.dna import DNA
from amplifyp.replication_origin import ReplicationOrigin
from amplifyp.replication_target import ReplicationConfig


def test_origin_init() -> None:
    """Test the initialization of a Origin object with invalid parameters."""
    with pytest.raises(ValueError):
        ReplicationOrigin(
            target="ATCG", primer="ATC", replication_config=ReplicationConfig()
        )


@dataclass(frozen=True, slots=True)
class OriginExample:
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


origin_examples: List[OriginExample] = []
origin_examples.append(
    OriginExample(
        ReplicationOrigin(
            target=DNA("CGACTGGGCAAAGGAAATCCTT").complement().reverse().sequence,
            primer="GCTGACCCNTTTCYYTTAGGCA",
            replication_config=ReplicationConfig(),
        ),
        0.9923605805958747,
        0.9118497898586322,
        0.8802629630681338,
    )
)

origin_examples.append(
    OriginExample(
        ReplicationOrigin(
            target=DNA("CGAGGGGGCAAAGGAAATCC").reverse().sequence,
            primer=DNA("CGACTGGGCAAAGGAAATCC").reverse().sequence,
            replication_config=ReplicationConfig(),
        ),
        0.9850746268656716,
        0.8914650537634409,
        0.8456746007863905,
    )
)


def test_origin_primability() -> None:
    """Test if origin primability is working correctly."""
    for ex in origin_examples:
        assert ex.primability == ex.origin.primability


def test_origin_stability() -> None:
    """Test if origin stability is working correctly."""
    for ex in origin_examples:
        assert ex.stability == ex.origin.stability


def test_origin_quality() -> None:
    """Test if origin quality is working correctly."""
    for ex in origin_examples:
        assert ex.quality == ex.origin.quality
