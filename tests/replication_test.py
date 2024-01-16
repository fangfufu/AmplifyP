# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNADirection, DNAType, Primer
from amplifyp.replication import Repliconf
from amplifyp.settings import Settings


def test_Repliconf_search() -> None:
    template = DNA("AAAAAATTTAAAAAAA")
    primer = Primer("TTT")
    repliconf = Repliconf(template, primer, settings=Settings)
