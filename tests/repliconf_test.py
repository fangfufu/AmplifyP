# -*- coding: utf-8 -*-
"""Tests related to repliconf."""

from amplifyp.dna import DNA, DNADirection, DNAType, Primer
from amplifyp.replication import Repliconf
from amplifyp.settings import DEFAULT_SETTINGS


def test_repliconf_linear_search() -> None:
    """A short test for the linear of repliconf."""
    # This is just for helping me to count the index
    #               0123456789ABCDEF
    template = DNA("ACCTCCTAGGAGGTTT")
    primer = Primer("CCT")
    repliconf = Repliconf(template, primer, DEFAULT_SETTINGS)
    assert repliconf.template_seq[DNADirection.FWD] == "---ACCTCCTAGGAGGTTT"
    assert repliconf.template_seq[DNADirection.REV] == "TGGAGGATCCTCCAAA---"
    repliconf.search()
    assert repliconf.origin_id[DNADirection.FWD] == [4, 7]
    assert repliconf.origin_id[DNADirection.REV] == [7, 10]
    assert repliconf.amplicon_start == [1, 4]
    assert repliconf.amplicon_end == [10, 13]


def test_repliconf_circular_search() -> None:
    """A short test for circular search of repliconf"""
    # This is just for helping me to count the index
    #               0123456789ABCDEF
    template = DNA("TGAAAAAGGAAAAACC", DNAType.CIRCULAR)
    primer = Primer("CCT")
    repliconf = Repliconf(template, primer, DEFAULT_SETTINGS)
    assert repliconf.template_seq[DNADirection.FWD] == "ACCTGAAAAAGGAAAAACC"
    assert repliconf.template_seq[DNADirection.REV] == "ACTTTTTCCTTTTTGGACT"
    repliconf.search()
    assert repliconf.origin_id[DNADirection.FWD] == [1]
    assert repliconf.origin_id[DNADirection.REV] == [6]
    assert repliconf.amplicon_start == [-2]
    assert repliconf.amplicon_end == [9]
