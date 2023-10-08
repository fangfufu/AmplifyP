# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNAType
from amplifyp.replisome import Replisome


def simple_replisome_test() -> None:
    """A simple replisome test."""
    target = DNA("GCTAGCGGAGTGTATACTGGCTTACT")
    primer = DNA("GCTA", DNAType.PRIMER)
    replisome = Replisome(target, primer)
    assert replisome.primability == 0
