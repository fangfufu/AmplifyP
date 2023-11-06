# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNADirection, DNAType, Primer
from amplifyp.replication import ReplicationConfig
from amplifyp.settings import Settings


def test_replication_config_linear() -> None:
    """Test the ReplicationConfig class with a linear target DNA."""
    target = DNA("ATCGATCG")
    primer = Primer("CGAT")
    min_overlap = 2
    config = Settings(min_overlap=min_overlap)

    rc = ReplicationConfig(target, primer, config)
    rc.search()
    assert primer.index[target, DNADirection.FWD] == [4]
    idx_start = primer.index[target, DNADirection.FWD][0]
    idx_end = idx_start + len(primer)
    assert (
        target.pad(len(primer) - min_overlap).sequence[idx_start:idx_end]
        == primer.sequence
    )
    assert not primer.index[target, DNADirection.REV]

    assert rc.primer == primer
    assert rc.primer_seq == "TAGC"
    assert rc.target_seq[DNADirection.FWD] == "GCTAGCTA--"
    assert rc.target_seq[DNADirection.REV] == "TAGCTAGC--"
    assert rc.range() == range(0, 7)
    assert rc.slice(2) == slice(2, 6)

    origin = rc.origin(DNADirection.FWD, 0)
    assert origin.primer == "TAGC"
    assert origin.target == "GCTA"

    origin = rc.origin(DNADirection.FWD, 6)
    assert origin.primer == "TAGC"
    assert origin.target == "TA--"

    origin = rc.origin(DNADirection.REV, 0)
    assert origin.primer == "TAGC"
    assert origin.target == "TAGC"

    origin = rc.origin(DNADirection.REV, 6)
    assert origin.primer == "TAGC"
    assert origin.target == "GC--"


def test_replication_config_circular() -> None:
    """Test the ReplicationConfig class with a circular target DNA."""
    target = DNA("ATCGATCG", DNAType.CIRCULAR)
    primer = Primer("CGAT")
    min_overlap = 1
    config = Settings(min_overlap=min_overlap)

    rc = ReplicationConfig(target, primer, config)
    rc.search()
    search_results = primer.index[target, DNADirection.FWD]
    assert search_results == [5, 1]
    for i in search_results:
        idx_start = i
        idx_end = idx_start + len(primer)
        assert (
            target.pad(len(primer) - min_overlap).sequence[idx_start:idx_end]
            == primer.sequence
        )
    assert not primer.index[target, DNADirection.REV]

    assert rc.primer == primer
    assert rc.primer_seq == "TAGC"
    assert rc.target_seq[DNADirection.FWD] == "GCTAGCTAGCT"
    assert rc.target_seq[DNADirection.REV] == "TAGCTAGCTAG"
    assert rc.range() == range(0, 8)
    assert rc.slice(2) == slice(2, 6)

    origin = rc.origin(DNADirection.FWD, 0)
    assert origin.primer == "TAGC"
    assert origin.target == "GCTA"

    origin = rc.origin(DNADirection.FWD, 6)
    assert origin.primer == "TAGC"
    assert origin.target == "TAGC"

    origin = rc.origin(DNADirection.REV, 0)
    assert origin.primer == "TAGC"
    assert origin.target == "TAGC"

    origin = rc.origin(DNADirection.REV, 6)
    assert origin.primer == "TAGC"
    assert origin.target == "GCTA"
