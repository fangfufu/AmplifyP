# -*- coding: utf-8 -*-
"""Simple tests for replisome.py."""

from amplifyp.dna import DNA, DNADirection, DNAType
from amplifyp.replication import ReplicationConfig
from amplifyp.settings import Settings


def test_replication_config_linear() -> None:
    """Test the ReplicationConfig class with a linear target DNA."""
    target = DNA("ATCGATCG")
    primer = DNA("CGAT")
    config = Settings(min_overlap=2)

    rc = ReplicationConfig(target, primer, config)

    assert rc.primer == primer
    assert rc.primer_sequence == "TAGC"
    assert rc.forward_sequence == "GCTAGCTA--"
    assert rc.reverse_sequence == "TAGCTAGC--"
    assert rc.range() == range(0, 7)
    assert rc.slice(2) == slice(2, 6)

    origin = rc.origin(DNADirection.FORWARD, 0)
    assert origin.primer == "TAGC"
    assert origin.target == "GCTA"

    origin = rc.origin(DNADirection.FORWARD, 6)
    assert origin.primer == "TAGC"
    assert origin.target == "TA--"

    origin = rc.origin(DNADirection.REVERSE, 0)
    assert origin.primer == "TAGC"
    assert origin.target == "TAGC"

    origin = rc.origin(DNADirection.REVERSE, 6)
    assert origin.primer == "TAGC"
    assert origin.target == "GC--"


def test_replication_config_circular() -> None:
    """Test the ReplicationConfig class with a circular target DNA."""
    target = DNA("ATCGATCG", DNAType.CIRCULAR)
    primer = DNA("CGAT")
    config = Settings(min_overlap=1)

    rc = ReplicationConfig(target, primer, config)

    assert rc.primer == primer
    assert rc.primer_sequence == "TAGC"
    assert rc.forward_sequence == "GCTAGCTAGCT"
    assert rc.reverse_sequence == "TAGCTAGCTAG"
    assert rc.range() == range(0, 8)
    assert rc.slice(2) == slice(2, 6)

    origin = rc.origin(DNADirection.FORWARD, 0)
    assert origin.primer == "TAGC"
    assert origin.target == "GCTA"

    origin = rc.origin(DNADirection.FORWARD, 6)
    assert origin.primer == "TAGC"
    assert origin.target == "TAGC"

    origin = rc.origin(DNADirection.REVERSE, 0)
    assert origin.primer == "TAGC"
    assert origin.target == "TAGC"

    origin = rc.origin(DNADirection.REVERSE, 6)
    assert origin.primer == "TAGC"
    assert origin.target == "GCTA"
