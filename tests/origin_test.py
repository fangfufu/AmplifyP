# -*- coding: utf-8 -*-
"""Tests related to replication origin."""
from dataclasses import dataclass
from typing import List
import pytest

from amplifyp.origin import ReplicationOrigin, Amplify4RevOrigin, Amplify4FwdOrigin
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

    origin: ReplicationOrigin | Amplify4RevOrigin | Amplify4FwdOrigin
    primability: float
    stability: float
    quality: float


origin_examples: List[ReplicationOriginExample] = []

# Origin example 1 from the documentation
origin_examples.append(
    ReplicationOriginExample(
        Amplify4RevOrigin(
            target="CGACTGGGCAAAGGAAATCCTT",
            primer="GCTGACCCNTTTCYYTTAGGCA",
        ),
        0.9923605805958747,  # Amplify4 reports 99%
        0.9293543192561425,  # Amplify4 reports 92%
        0.9021436248150215,  # Amplify4 reports 0.8875
    )
)

# Origin example 2 from the documentation
origin_examples.append(
    ReplicationOriginExample(
        Amplify4FwdOrigin(
            target="CGAGGGGGCAAAGGAAATCC",
            primer="CGACTGGGCAAAGGAAATCC",
        ),
        0.9850746268656716,  # Amplify4 reports 98%
        0.8911290322580645,  # Amplify4 reports 89%
        0.8452545739046701,  # Amplify4 reports 0.8375
    )
)

# Reverse Origin for test_verification_F1423_R20049, Primer 20049
origin_examples.append(
    ReplicationOriginExample(
        Amplify4RevOrigin(
            target="ATCCTCTTTTTATGTCCACATCTG",
            primer="TAACAGAAAAATACAGGTGTAGAC",
        ),
        0.855072463768116,  # Amplify4 reports 85%
        0.9005376344086021,  # Amplify4 reports 90%
        0.6945126227208976,  # Amplify4 reports 0.6875
    )
)

# Forward Origin for test_verification_F1423_R20049, Primer 1423
origin_examples.append(
    ReplicationOriginExample(
        Amplify4FwdOrigin(
            target="CGCAAAGCTTGTCGGC",
            primer="GCCAAAGTGTGTTGGC",
        ),
        0.8076923076923077,  # Amplify4 reports 80%
        0.6653225806451613,  # Amplify4 reports 66%
        0.34126861042183604,  # Amplify4 reports 0.325
    )
)

origin_examples.append(
    ReplicationOriginExample(
        Amplify4FwdOrigin(
            target="GAAAAAATTAACAAAAAAAACAAT",
            primer="CAGATGTGGACATAAAAAGACAAT",
        ),
        0.8188405797101449,  # Amplify4 reports 81%
        0.4838709677419355,  # Amplify4 reports 48%
        0.12838943431510044,  # Amplify4 reports 0.1125
    )
)

origin_examples.append(
    ReplicationOriginExample(
        Amplify4FwdOrigin(
            target="CCTTTTT",
            primer="TTTTTTT",
        ),
        0.8229166666666666,  # Amplify4 reports 82%
        0.7142857142857143,  # Amplify4 reports 71%
        0.42150297619047605,  # Amplify4 reports 0.4125
    )
)

origin_examples.append(
    ReplicationOriginExample(
        Amplify4FwdOrigin(
            target="CTTTTTT",
            primer="TTTTTTT",
        ),
        0.9166666666666666,  # Amplify4 reports 91%
        0.8571428571428571,  # Amplify4 reports 85%
        0.7172619047619045,  # Amplify4 reports 0.7
    )
)

origin_examples.append(
    ReplicationOriginExample(
        Amplify4FwdOrigin(
            target="TTTTTTT",
            primer="TTTTTTT",
        ),
        1.0,  # Amplify4 reports 100%
        1.0,  # Amplify4 reports 100%
        1.0,  # Amplify4 reports 1.0
    )
)

origin_examples.append(
    ReplicationOriginExample(
        Amplify4FwdOrigin(
            target="CCTTTTTT",
            primer="TTTTTTTT",
        ),
        0.8543689320388349,  # Amplify4 reports 85%
        0.75,  # Amplify4 reports 75%
        0.5054611650485437,  # Amplify4 reports 0.5
    )
)

origin_examples.append(
    ReplicationOriginExample(
        Amplify4FwdOrigin(
            target="CTTTTTTT",
            primer="TTTTTTTT",
        ),
        0.9320388349514563,  # Amplify4 reports 93%
        0.875,  # Amplify4 reports 87%
        0.7587985436893204,  # Amplify4 reports 0.75
    )
)

origin_examples.append(
    ReplicationOriginExample(
        Amplify4FwdOrigin(
            target="TTTTTTTT",
            primer="TTTTTTTT",
        ),
        1.0,  # Amplify4 reports 100%
        1.0,  # Amplify4 reports 100%
        1.0,  # Amplify4 reports 1.0
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
        assert ex.stability == ex.origin.stability()


def test_origin_quality() -> None:
    """Test if origin quality is working correctly."""
    print("")
    for ex in origin_examples:
        assert ex.quality == ex.origin.quality()
