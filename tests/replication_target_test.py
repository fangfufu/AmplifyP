# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from dataclasses import dataclass
from typing import List

from amplifyp.dna import DNA, DNAType


@dataclass(frozen=True, slots=True)
class ReplicationTargetExample:
    """A class representing test cases of replisome.

    Attributes:
        target (DNA): The target DNA sequence.
        primer (DNA): The primer DNA sequence.
    """

    target: DNA
    primer: DNA
    origin_index: int


replisome_examples: List[ReplicationTargetExample] = []
replisome_examples.append(
    ReplicationTargetExample(
        target=DNA("GTGTACTCAGTTCCATAAACGAGCTATTAGATATGAGGTCCGTAGATTGAAAAGGGTGA"),
        primer=DNA("TTCCATA", DNAType.PRIMER),
        origin_index=10,
    )
)

replisome_examples.append(
    ReplicationTargetExample(
        target=DNA("ATTGTGCGATCCCTGCACCTCAGCTAAGGTAGCTACCAATATTTAGTTTCTAAGCCTTGC"),
        primer=DNA("GCGATC", DNAType.PRIMER),
        origin_index=5,
    )
)
