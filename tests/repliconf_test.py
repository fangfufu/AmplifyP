# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNADirection, Primer
from amplifyp.replication import Repliconf
from amplifyp.settings import Settings


def test_repliconf_search_short() -> None:
    """A short test for the search feature of repliconf."""
    # This is just for helping me to count the index
    #               0123456789ABCDEF
    template = DNA("ACCTCCTAGGAGGTTT")
    primer = Primer("CCT")
    repliconf = Repliconf(template, primer, Settings)
    assert repliconf.template_seq[DNADirection.FWD] == "---ACCTCCTAGGAGGTTT"
    assert repliconf.template_seq[DNADirection.REV] == "TGGAGGATCCTCCAAA---"
    repliconf.search()
    assert repliconf.origin_id[DNADirection.FWD] == [4, 7]
    assert repliconf.origin_id[DNADirection.REV] == [7, 10]
    assert repliconf.amplicon_start == [1, 4]
    assert repliconf.amplicon_end == [10, 13]
