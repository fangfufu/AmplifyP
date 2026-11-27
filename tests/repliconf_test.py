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

from amplifyp.dna import DNA, DNADirection, DNAType, Nucleotides, Primer
from amplifyp.origin import ReplicationOrigin
from amplifyp.repliconf import DirIdx, Repliconf
from amplifyp.settings import (
    GLOBAL_REPLICATION_SETTINGS,
    BasePairWeightsTbl,
    LengthWiseWeightTbl,
    ReplicationSettings,
)


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


def test_settings_wrapper_equality() -> None:
    """Test _SettingsWrapper equality with non-_SettingsWrapper objects."""
    from amplifyp.repliconf import _SettingsWrapper
    from amplifyp.settings import ReplicationSettings

    wrapper = _SettingsWrapper(ReplicationSettings())
    assert wrapper != "not_a_wrapper"
    assert wrapper != 123


def test_repliconf_origin_from_db() -> None:
    """Test Repliconf.origin_from_db lookup by direction and search index."""
    test_repliconf.search()
    orig_fwd = test_repliconf.origin_from_db(DNADirection.FWD, 0)
    assert orig_fwd is not None
    orig_rev = test_repliconf.origin_from_db(DNADirection.REV, 0)
    assert orig_rev is not None


def test_repliconf_format_context_alignment() -> None:
    """Test Repliconf.format_context_alignment for FWD and REV origins."""
    template = DNA("ACCTCCTAGGAGGTTT", DNAType.LINEAR)
    primer = Primer("CCT", name="MyPrimer")
    conf = Repliconf(template, primer, GLOBAL_REPLICATION_SETTINGS)
    conf.search()

    fwd_var = conf.origin_db.fwd[0]
    fwd_align = conf.format_context_alignment(
        fwd_var, improved_visualisation=True, separator="-"
    )
    assert "5'-CCT-3'" in fwd_align.primer_line
    assert fwd_align.complement_line is not None
    assert "3'-" in fwd_align.complement_line
    assert "5'-" in fwd_align.template_line
    assert len(fwd_align.top_block.splitlines()) == 3
    assert len(fwd_align.bottom_block.splitlines()) == 3
    assert len(fwd_align.to_text().splitlines()) == 7

    rev_var = conf.origin_db.rev[0]
    rev_align = conf.format_context_alignment(
        rev_var, improved_visualisation=True, separator=" "
    )
    assert "MyPrimer (Reverse)" in rev_align.primer_line
    assert "3' TCC 5'" in rev_align.primer_line
    assert rev_align.complement_line is None
    assert len(rev_align.bottom_block.splitlines()) == 2

    empty_conf = Repliconf(DNA(""), Primer("CCT"))
    with pytest.raises(ValueError, match="empty template or primer"):
        empty_conf.format_context_alignment(DirIdx(DNADirection.FWD, 0))

    # Variable context_bases alignment keeps primer and context columns aligned
    ctx10 = conf.format_context_alignment(fwd_var, context_bases=10)
    assert ctx10.primer_line.index("CCT") == ctx10.template_line.index("CCT")
    assert ctx10.pipe_line.index("|") == ctx10.primer_line.index("CCT")
    assert ctx10.arrow_line.index("V") == ctx10.primer_line.index("CCT")

    # 1-bp primer renders single marker and validates invalid arguments
    one_bp_conf = Repliconf(DNA("ACGT"), Primer("C", name="P1"))
    one_bp_conf.search()
    one_bp_align = one_bp_conf.format_context_alignment(
        one_bp_conf.origin_db.fwd[0], context_bases=0
    )
    assert one_bp_align.pipe_line.strip() == "|"
    assert one_bp_align.arrow_line.strip() == "V"

    with pytest.raises(ValueError, match="separator"):
        conf.format_context_alignment(fwd_var, separator="Dash ('-')")
    with pytest.raises(ValueError, match="context_bases"):
        conf.format_context_alignment(fwd_var, context_bases=-1)

    # Linear 5' overhang gap-pads 5' end rather than wrapping to 3' end
    linear_tpl = DNA("AAAAACCCCCTTTTTGGGGG", DNAType.LINEAR)
    overhang_primer = Primer("TTAAAA", name="Overhang")
    overhang_conf = Repliconf(linear_tpl, overhang_primer)
    overhang_align = overhang_conf.format_context_alignment(
        DirIdx(DNADirection.FWD, 4), context_bases=2
    )
    assert "5' ----AAAAAC 3'" in overhang_align.template_line
    assert "-2" in overhang_align.coordinates_line
    assert "4" in overhang_align.coordinates_line

    # Circular template wraps context slices and coordinates cleanly
    circ_tpl = DNA("ACGTACGTACGT", DNAType.CIRCULAR)
    circ_conf = Repliconf(circ_tpl, Primer("GTAC", name="Circ"))
    circ_align = circ_conf.format_context_alignment(
        DirIdx(DNADirection.FWD, 6), context_bases=4
    )
    assert "5' GTACGTACGTAC 3'" in circ_align.template_line

    # Short 3-bp primer at 3-digit coordinates with context_bases=0 does not
    # overwrite start/end labels or clip trailing digits
    long_tpl = DNA("A" * 120, DNAType.LINEAR)
    short_conf = Repliconf(long_tpl, Primer("AAA", name="Short"))
    short_align = short_conf.format_context_alignment(
        DirIdx(DNADirection.FWD, 102), context_bases=0
    )
    assert "100 102" in short_align.coordinates_line


def test_repliconf_pruning_fallback_and_lazy_properties() -> None:
    """Test pruning fallback paths, degenerate bases, and lazy reverse comp."""
    # Long 37-mer primer where a 3'-end mismatch still scores 121/151 > 0.8
    # (good_acgt_count > 2 -> triggers non-masked full search_range scan).
    long_primer = Primer("A" * 37)
    tpl = DNA("TTTT" + ("A" * 36) + "CTTTT", DNAType.LINEAR)
    conf = Repliconf(tpl, long_primer)
    conf.search()
    mismatch_idx = len(long_primer) + 4  # 37 gap bases + "TTTT" -> "A"*36 + "C"
    mismatch_origins = [
        o for o in conf.origin_db.fwd if o.index == mismatch_idx
    ]
    assert len(mismatch_origins) == 1
    assert mismatch_origins[0].primability == pytest.approx(121 / 151)
    assert conf.rev_comp_primer == long_primer.reverse_complement()

    # Non-positive primability cutoff / lower cutoff exercises rem_max & low cut
    permissive = ReplicationSettings(
        primability_cutoff=-0.2,
        stability_cutoff=0.1,
        match_weight=LengthWiseWeightTbl(1.0),
        run_weights=LengthWiseWeightTbl(100.0),
    )
    perm_conf = Repliconf(
        DNA("acgtACGT", DNAType.LINEAR), Primer("ncgt"), permissive
    )
    perm_conf.search()
    assert len(perm_conf.origin_db.fwd) > 0

    # Empty and shorter-than-primer circular templates must not crash
    empty_conf = Repliconf(DNA("", DNAType.LINEAR), Primer("ACGT"))
    empty_conf.search()
    assert len(empty_conf.origin_db.fwd) == 0

    tiny_circ = Repliconf(DNA("ACGT", DNAType.CIRCULAR), Primer("ACGTACGTACGT"))
    tiny_circ.search()
    assert len(tiny_circ.origin_db.fwd) == 0

    # search() does not eagerly materialise template_seq
    lazy_conf = Repliconf(DNA("ACGTACGT", DNAType.LINEAR), Primer("ACGT"))
    lazy_conf.search()
    assert "template_seq" not in lazy_conf.__dict__
    assert len(lazy_conf.origin_db.fwd) >= 1
    assert lazy_conf.template_seq[DNADirection.FWD] == "----ACGTACGT"


def test_repliconf_search_pruning_parity_with_origin_bruteforce() -> None:
    """Verify masked/pruned search matches brute-force ReplicationOrigin."""
    # Custom weights where only 'N' passes at 3' end (good_acgt_count == 0
    # while has_any_good_k0 is True).
    only_n_matrix = [[0.0, 0.0, 0.0, 0.0, 100.0] for _ in Nucleotides.PRIMER]
    only_n_settings = ReplicationSettings(
        base_pair_scores=BasePairWeightsTbl(
            Nucleotides.PRIMER, Nucleotides.TEMPLATE, only_n_matrix
        ),
        primability_cutoff=0.8,
        stability_cutoff=0.4,
    )

    cases: list[tuple[DNA, Primer, ReplicationSettings]] = [
        (
            DNA("GGACGTACGTACGTNNACGTACGTACGTCC", DNAType.LINEAR),
            Primer("ACGTACGT"),
            GLOBAL_REPLICATION_SETTINGS,
        ),
        (
            DNA("ACGTACGTNNTTGGCCAA", DNAType.CIRCULAR),
            Primer("ACGTMRWS"),
            GLOBAL_REPLICATION_SETTINGS,
        ),
        (
            DNA("ACGTTTTNACGT", DNAType.LINEAR),
            Primer("ACGT"),
            ReplicationSettings(primability_cutoff=0.5, stability_cutoff=0.2),
        ),
        (
            DNA("ACGTNNNNCGTA", DNAType.LINEAR),
            Primer("ACGT"),
            only_n_settings,
        ),
    ]

    for tpl, primer, settings in cases:
        conf = Repliconf(tpl, primer, settings)
        conf.search()
        for direction in (DNADirection.FWD, DNADirection.REV):
            expected: list[tuple[int, float, float]] = []
            for i in conf.range():
                orig = conf.origin(DirIdx(direction, i))
                ref = ReplicationOrigin(orig.target, orig.primer, settings)
                if (
                    ref.primability > settings.primability_cutoff
                    and ref.stability > settings.stability_cutoff
                ):
                    expected.append((i, ref.primability, ref.stability))

            actual = (
                conf.origin_db.fwd
                if direction == DNADirection.FWD
                else conf.origin_db.rev
            )
            assert [o.index for o in actual] == [idx for idx, _, _ in expected]
            for obs, (_, exp_p, exp_s) in zip(actual, expected, strict=True):
                assert obs.primability == pytest.approx(exp_p)
                assert obs.stability == pytest.approx(exp_s)


def test_repliconf_search_reflects_in_place_table_mutations() -> None:
    """Verify in-place edits to weight tables invalidate cached scoring."""
    # 1. In-place mutation of run_weights
    rw_settings = ReplicationSettings(
        primability_cutoff=0.4,
        stability_cutoff=0.4,
        match_weight=LengthWiseWeightTbl(1.0),
        run_weights=LengthWiseWeightTbl(1.0),
    )
    rw_conf = Repliconf(DNA("AC", DNAType.LINEAR), Primer("AT"), rw_settings)
    rw_conf.search()
    assert any(o.index == 2 for o in rw_conf.origin_db.fwd)

    rw_settings.run_weights[0] = 0.0
    rw_conf.search()
    assert not any(o.index == 2 for o in rw_conf.origin_db.fwd)

    # 2. In-place mutation of match_weight
    mw_settings = ReplicationSettings(
        primability_cutoff=0.4,
        stability_cutoff=0.4,
        match_weight=LengthWiseWeightTbl(1.0),
        run_weights=LengthWiseWeightTbl(1.0),
    )
    mw_conf = Repliconf(DNA("AC", DNAType.LINEAR), Primer("AT"), mw_settings)
    mw_conf.search()
    assert any(o.index == 2 for o in mw_conf.origin_db.fwd)

    mw_settings.match_weight[1] = 0.0
    mw_conf.search()
    assert not any(o.index == 2 for o in mw_conf.origin_db.fwd)

    # 3. In-place mutation of base_pair_scores
    bp_settings = ReplicationSettings(
        primability_cutoff=0.4,
        stability_cutoff=0.4,
        match_weight=LengthWiseWeightTbl(1.0),
        run_weights=LengthWiseWeightTbl(1.0),
    )
    bp_conf = Repliconf(DNA("AC", DNAType.LINEAR), Primer("AT"), bp_settings)
    bp_conf.search()
    assert any(o.index == 2 for o in bp_conf.origin_db.fwd)

    bp_settings.base_pair_scores["A", "A"] = 0.0
    bp_conf.search()
    assert not any(o.index == 2 for o in bp_conf.origin_db.fwd)


def test_repliconf_template_seq_access_does_not_alter_search() -> None:
    """Verify materialising template_seq keeps non-ASCII search consistent."""
    for dna_type in (DNAType.LINEAR, DNAType.CIRCULAR):
        tpl = DNA("A\u017f", dna_type)
        primer = Primer("A")

        cold_conf = Repliconf(tpl, primer)
        cold_conf.search()

        warm_conf = Repliconf(tpl, primer)
        _ = warm_conf.template_seq
        warm_conf.search()

        assert warm_conf.origin_db.fwd == cold_conf.origin_db.fwd
        assert warm_conf.origin_db.rev == cold_conf.origin_db.rev
        assert len(warm_conf.origin_db.fwd) == 1

    # Preserves explicit cached template_seq overrides with uppercase encoding
    override_conf = Repliconf(DNA("TT", DNAType.LINEAR), Primer("A"))
    override_conf.template_seq = {
        DNADirection.FWD: "-a\u017f",
        DNADirection.REV: "t\u017f-",
    }
    override_conf.search()
    assert [o.index for o in override_conf.origin_db.fwd] == [1]

    # Circular non-ASCII padding uses raw coordinates before uppercase expansion
    unit_settings = ReplicationSettings(
        primability_cutoff=0.4,
        stability_cutoff=0.4,
        match_weight=LengthWiseWeightTbl(1.0),
        run_weights=LengthWiseWeightTbl(1.0),
    )
    for seq_str, primer_str, exp_fwd_bytes, exp_rev_bytes in (
        ("Aß", "A", b"SSASS", b"TSST"),
        ("ßA", "A", b"ASSA", b"SSTSS"),
        ("A\ufb05", "T", b"STAST", b"TSTT"),
    ):
        circ_tpl = DNA(seq_str, DNAType.CIRCULAR)
        circ_primer = Primer(primer_str)
        cold_circ = Repliconf(circ_tpl, circ_primer, unit_settings)
        assert cold_circ._get_template_bytes(DNADirection.FWD) == exp_fwd_bytes
        assert cold_circ._get_template_bytes(DNADirection.REV) == exp_rev_bytes
        cold_circ.search()

        warm_circ = Repliconf(circ_tpl, circ_primer, unit_settings)
        _ = warm_circ.template_seq
        warm_circ.search()

        assert cold_circ.origin_db.fwd == warm_circ.origin_db.fwd
        assert cold_circ.origin_db.rev == warm_circ.origin_db.rev
        if seq_str == "Aß":
            assert [o.index for o in cold_circ.origin_db.fwd] == [2]

    # Circular template with empty primer exercises zero-length padding path
    zero_pad_circ = Repliconf(DNA("ACGT", DNAType.CIRCULAR), Primer(""))
    zero_pad_circ.search()
    assert zero_pad_circ.searched
    assert len(zero_pad_circ.origin_db.fwd) == 0


def test_repliconf_rejects_length_changing_primer_normalisation() -> None:
    """Reject primers whose uppercase form changes coordinate length."""
    permissive = ReplicationSettings(
        primability_cutoff=-0.2,
        stability_cutoff=-0.2,
        match_weight=LengthWiseWeightTbl(1.0),
        run_weights=LengthWiseWeightTbl(1.0),
    )
    tpl = DNA("GGCC", DNAType.LINEAR)
    expanding_primer = Primer("ß")

    with pytest.raises(ValueError, match="uppercase normalisation"):
        Repliconf(tpl, expanding_primer, permissive)

    # Also reject if primer is swapped after construction before search()
    conf = Repliconf(tpl, Primer("G"), permissive)
    conf.primer = expanding_primer
    with pytest.raises(ValueError, match="uppercase normalisation"):
        conf.search()
    assert not conf.searched
    assert len(conf.origin_db.fwd) == 0
    assert len(conf.origin_db.rev) == 0


def test_repliconf_zero_denom_and_impossible_cutoff_paths() -> None:
    """Cover zero prim_denom and all-zero k0_mask early-skip branches."""
    # Zero match weights -> prim_denom == 0 -> unmasked fallback with 0 hits
    zero_settings = ReplicationSettings(
        match_weight=LengthWiseWeightTbl(0.0),
    )
    zero_conf = Repliconf(DNA("ACGTACGT"), Primer("ACGT"), zero_settings)
    zero_conf.search()
    assert zero_conf.searched
    assert len(zero_conf.origin_db.fwd) == 0
    assert len(zero_conf.origin_db.rev) == 0

    # Cutoff > 1.0 -> no k0 base can reach prim_req -> skips both strands
    impossible_settings = ReplicationSettings(primability_cutoff=1.5)
    impossible_conf = Repliconf(
        DNA("ACGTACGT"), Primer("ACGT"), impossible_settings
    )
    impossible_conf.search()
    assert impossible_conf.searched
    assert len(impossible_conf.origin_db.fwd) == 0
    assert len(impossible_conf.origin_db.rev) == 0
