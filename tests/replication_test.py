# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNADirection, Primer
from amplifyp.replication import Repliconf
from amplifyp.settings import Settings


def test_repliconf_search_short() -> None:
    """A short test for the search feature of repliconf."""
    # This is just for helping me to count the index
    #               0123456789ABCDEF
    template = DNA("AACCTCCTAGGAGGTT")
    primer = Primer("CCT")
    repliconf = Repliconf(template, primer, Settings)
    repliconf.search()
    assert repliconf.origin_id[DNADirection.FWD] == [5, 8]
    assert repliconf.origin_id[DNADirection.REV] == [8, 11]
    assert repliconf.amplicon_start == [2, 5]
    assert repliconf.amplicon_end == [11, 14]
