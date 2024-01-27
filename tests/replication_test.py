# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNADirection, DNAType, Primer
from amplifyp.replication import Repliconf
from amplifyp.settings import Settings


def test_repliconf_search_short() -> None:
    """A short test for the search feature of repliconf."""
    # This is just for helping me to count the index
    #               0123456789ABCDEF
    template = DNA("AACCTCCTAGGAGGTT")
    primer = Primer("CCT")
    repliconf = Repliconf(template, primer, settings=Settings)
    repliconf.search()
    print(repliconf.amplicon_start)
    print(repliconf.amplicon_end)
