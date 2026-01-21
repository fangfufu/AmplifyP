"""Tests related to repliconf."""

from amplifyp.dna import DNA, DNADirection, DNAType, Primer
from amplifyp.repliconf import Repliconf
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
    assert repliconf.origin_db.fwd == [4, 7]
    assert repliconf.origin_db.rev == [7, 10]


def test_repliconf_circular_search() -> None:
    """A short test for circular search of repliconf."""
    # This is just for helping me to count the index
    #               0123456789ABCDEF
    template = DNA("TGAAAAAGGAAAAACC", DNAType.CIRCULAR)
    primer = Primer("CCT")
    repliconf = Repliconf(template, primer, DEFAULT_SETTINGS)
    assert repliconf.template_seq[DNADirection.FWD] == "ACCTGAAAAAGGAAAAACC"
    assert repliconf.template_seq[DNADirection.REV] == "ACTTTTTCCTTTTTGGACT"
    repliconf.search()
    assert repliconf.origin_db.fwd == [1]
    assert repliconf.origin_db.rev == [6]


test_repliconf = Repliconf(
    DNA("TGAAAAAGGAAAAACC", DNAType.CIRCULAR), Primer("CCT"), DEFAULT_SETTINGS
)


def test_repliconf_comparison() -> None:
    """Test repliconf comparison."""
    a = test_repliconf
    assert a == test_repliconf
    assert test_repliconf != ""


def test_repliconf_str() -> None:
    """Test repliconf string representation."""
    assert str(test_repliconf) == (
        "ReplicationConfig: Primer: DNA: CCT, "
        + "PRIMER, FWD, Target: DNA: TGAAAAAGGAAAAACC, CIRCULAR, FWD"
    )


def test_repliconf_hash() -> None:
    """Test that the Repliconf hash function generates a hash."""
    assert hash(test_repliconf)
