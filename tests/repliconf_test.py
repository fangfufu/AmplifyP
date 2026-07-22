# Copyright (C) 2026 AmplifyP Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Tests related to repliconf."""

import pytest

from amplifyp.dna import DNA, DNADirection, DNAType, Primer
from amplifyp.repliconf import DirIdx, Repliconf
from amplifyp.settings import GLOBAL_REPLICATION_SETTINGS


def test_repliconf_linear_search() -> None:
    """A short test for the linear of repliconf."""
    # This is just for helping me to count the index
    #               0123456789ABCDEF
    template = DNA("ACCTCCTAGGAGGTTT")
    primer = Primer("CCT")
    repliconf = Repliconf(template, primer, GLOBAL_REPLICATION_SETTINGS)
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
    repliconf = Repliconf(template, primer, GLOBAL_REPLICATION_SETTINGS)
    assert repliconf.template_seq[DNADirection.FWD] == "ACCTGAAAAAGGAAAAACC"
    assert repliconf.template_seq[DNADirection.REV] == "ACTTTTTCCTTTTTGGACT"
    repliconf.search()
    assert repliconf.origin_db.fwd == [1]
    assert repliconf.origin_db.rev == [6]


test_repliconf = Repliconf(
    DNA("TGAAAAAGGAAAAACC", DNAType.CIRCULAR),
    Primer("CCT"),
    GLOBAL_REPLICATION_SETTINGS,
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


def test_dir_idx_methods() -> None:
    """Test methods of DirIdx."""
    d1 = DirIdx(DNADirection.FWD, 10)
    d2 = DirIdx(DNADirection.FWD, 5)  # Same direction
    d3 = DirIdx(DNADirection.REV, 10)  # Different direction

    # __int__ and __index__
    assert int(d1) == 10
    assert list(range(20))[d1] == 10

    # __add__
    assert (d1 + 5).index == 15
    assert (d1 + d2).index == 15
    with pytest.raises(TypeError):
        _ = d1 + "a"

    # __sub__
    assert (d1 - 5).index == 5
    assert (d1 - d2).index == 5
    with pytest.raises(TypeError):
        _ = d1 - "a"

    # __eq__
    assert d1 == 10
    assert d1 == DirIdx(DNADirection.FWD, 10)
    assert d1 != d2
    assert d1 != d3  # diff direction
    assert d1 != "invalid"

    # Comparisons
    # lt
    assert d2 < d1
    assert d2 < 10
    with pytest.raises(TypeError):
        _ = d1 < "a"

    # gt
    assert d1 > d2
    assert d1 > 5
    with pytest.raises(TypeError):
        _ = d1 > "a"

    # le
    assert d2 <= d1
    assert d1 <= d1  # reflexivity
    assert d2 <= 10
    with pytest.raises(TypeError):
        _ = d1 <= "a"

    # ge
    assert d1 >= d2
    assert d1 >= d1  # reflexivity
    assert d1 >= 5
    with pytest.raises(TypeError):
        _ = d1 >= "a"

    # __str__
    assert str(d1) == "10"


def test_repliconf_settings_equality() -> None:
    """Test that Repliconf equality takes ReplicationSettings into account."""
    from amplifyp.settings import ReplicationSettings

    s1 = ReplicationSettings(primability_cutoff=0.5)
    s2 = ReplicationSettings(primability_cutoff=0.8)
    template = DNA("TGAAAAAGGAAAAACC", DNAType.CIRCULAR)
    primer = Primer("CCT")
    r1 = Repliconf(template, primer, s1)
    r2 = Repliconf(template, primer, s2)
    assert r1 != r2
