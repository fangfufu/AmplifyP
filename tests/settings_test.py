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

"""Tests related to settings.py."""

import pytest

from amplifyp.errors import ColumnLengthMismatchError, RowLengthMismatchError
from amplifyp.settings import (
    BasePairWeightsTbl,
    LengthWiseWeightTbl,
    PrimerDimerSettings,
    ReplicationSettings,
    TMSettings,
)


def test_length_wise_weight_tbl() -> None:
    """Test the LengthWiseWeightTbl class."""
    weight_tbl = LengthWiseWeightTbl(0.5, {0: 0.6, 1: 0.7})
    assert weight_tbl[0] == pytest.approx(0.6)
    assert weight_tbl[1] == pytest.approx(0.7)
    weight_tbl[2] = 0.8
    assert weight_tbl[2] == pytest.approx(0.8)
    assert weight_tbl[100] == pytest.approx(0.5)
    assert weight_tbl.default_weight == pytest.approx(0.5)
    assert weight_tbl.overrides == {0: 0.6, 1: 0.7, 2: 0.8}


def test_length_wise_weight_tbl_eq_and_repr() -> None:
    """Test equality and string representation of LengthWiseWeightTbl."""
    tbl1 = LengthWiseWeightTbl(0.5, {0: 0.6})
    tbl2 = LengthWiseWeightTbl(0.5, {0: 0.6})
    tbl3 = LengthWiseWeightTbl(0.5, {0: 0.7})

    assert tbl1 == tbl2
    assert tbl1 != tbl3
    assert tbl1 != "not_a_tbl"
    assert (
        repr(tbl1)
        == "LengthWiseWeightTbl(default_weight=0.5, overrides={0: 0.6})"
    )


def test_empty_length_wise_weight_tbl() -> None:
    """Test the initialisation of an empty LengthWiseWeightTbl."""
    a = LengthWiseWeightTbl()
    assert a[0] == 0


pairwise_weights = [
    [1.0, 0.5, 0.2, 0.1],
    [0.5, 1.0, 0.1, 0.2],
    [0.2, 0.1, 0.2, 0.5],
    [0.1, 0.2, 0.5, 0.2],
]


def test_invalid_tbl_generation() -> None:
    """Test invalid BasePairWeightsTbl creation."""
    with pytest.raises(RowLengthMismatchError):
        BasePairWeightsTbl("AB", "ABCD", pairwise_weights)

    with pytest.raises(ColumnLengthMismatchError):
        BasePairWeightsTbl("ABCD", "AB", pairwise_weights)


def test_base_pair_weights_tbl() -> None:
    """Test function for the BasePairWeightsTbl class."""
    # Define the nucleotides and their pairwise weights
    nucleotides = "ACGT-"

    # Create a BasePairWeightsTbl instance
    npwt = BasePairWeightsTbl(nucleotides, nucleotides, pairwise_weights)

    # Test the row and column properties
    assert npwt.row() == "ACGT"
    assert npwt.column() == "ACGT"

    # Test the __setitem__ method
    npwt[("G", "G")] = 1.0
    npwt[("T", "T")] = 1.0

    # Test the __getitem__ method
    assert npwt[("A", "C")] == pytest.approx(0.5)
    assert npwt[("G", "T")] == pytest.approx(0.5)
    assert npwt[("C", "A")] == pytest.approx(0.5)
    assert npwt[("T", "G")] == pytest.approx(0.5)
    assert npwt[("A", "A")] == pytest.approx(1.0)
    assert npwt[("C", "C")] == pytest.approx(1.0)
    assert npwt[("G", "G")] == pytest.approx(1.0)
    assert npwt[("T", "T")] == pytest.approx(1.0)

    # Test the support for gap symbol
    assert npwt["-", "A"] == 0

    assert (
        str(npwt)
        == "{('A', 'A'): 1.0, ('A', 'C'): 0.5, ('A', 'G'): 0.2, \
('A', 'T'): 0.1, ('A', '-'): 0, ('C', 'A'): 0.5, ('C', 'C'): 1.0, ('C', 'G'): \
0.1, ('C', 'T'): 0.2, ('C', '-'): 0, ('G', 'A'): 0.2, ('G', 'C'): 0.1, \
('G', 'G'): 1.0, ('G', 'T'): 0.5, ('G', '-'): 0, ('T', 'A'): 0.1, ('T', 'C'): \
0.2, ('T', 'G'): 0.5, ('T', 'T'): 1.0, ('T', '-'): 0, ('-', 'A'): 0, \
('-', 'C'): 0, ('-', 'G'): 0, ('-', 'T'): 0, ('-', '-'): 0}"
    )

    assert len(npwt) == 16


def test_base_pair_weights_tbl_case_sensitivity() -> None:
    """Test that BasePairWeightsTbl handles case sensitivity correctly."""
    row = "acgt-"
    col = "acgt-"
    weights = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]

    # Initialise with lowercase
    tbl = BasePairWeightsTbl(row, col, weights)

    # Check if row() and column() return uppercase (accessing row and col)
    assert tbl.row() == "ACGT"
    assert tbl.column() == "ACGT"

    # Check if __getitem__ works with uppercase keys
    assert tbl["A", "A"] == pytest.approx(1.0)

    # Check if __getitem__ works with lowercase keys
    assert tbl["a", "a"] == pytest.approx(1.0)

    # Check mixed case
    assert tbl["A", "a"] == pytest.approx(1.0)


def test_base_pair_weights_tbl_row_max_update() -> None:
    """Test that row_max is updated when weights are modified."""
    row = "ACGT-"
    col = "ACGT-"
    tbl = BasePairWeightsTbl(row, col, pairwise_weights)

    assert tbl.row_max("A") == pytest.approx(1.0)
    tbl["A", "T"] = 2.5
    assert tbl.row_max("A") == pytest.approx(2.5)


def test_base_pair_weights_tbl_without_gap_suffix() -> None:
    """Test BasePairWeightsTbl when row/col do not end with a gap symbol."""
    row = "ACGT"
    col = "ACGT"
    weights = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    tbl = BasePairWeightsTbl(row, col, weights)
    assert tbl.row() == "ACGT"
    assert tbl.column() == "ACGT"


def test_base_pair_weights_tbl_eq_and_repr() -> None:
    """Test equality and string representation of BasePairWeightsTbl."""
    row = "ACGT-"
    col = "ACGT-"
    tbl1 = BasePairWeightsTbl(row, col, pairwise_weights)
    tbl2 = BasePairWeightsTbl(row, col, pairwise_weights)

    assert tbl1 == tbl2
    assert repr(tbl1) == "BasePairWeightsTbl(row='ACGT', col='ACGT')"


def test_base_pair_weights_tbl_invalid_chars() -> None:
    """Test BasePairWeightsTbl with invalid characters (code >= 128)."""
    weird_char = "\u00c0"
    row = "A" + weird_char
    col = "A"
    weights = [[1.0], [2.0]]

    tbl = BasePairWeightsTbl(row, col, weights)

    # Should work via fallback or expanded map
    assert tbl[weird_char, "A"] == pytest.approx(2.0)

    # Setter should also work
    tbl[weird_char, "A"] = 3.0
    assert tbl[weird_char, "A"] == pytest.approx(3.0)

    # Also test completely invalid lookup
    with pytest.raises(KeyError):
        _ = tbl["Z", "Z"]


def test_tm_settings_dntp_conc() -> None:
    """Test dntp_conc field in TMSettings."""
    tm = TMSettings(dntp_conc=2.5)
    assert tm.dntp_conc == pytest.approx(2.5)
    tm.dntp_conc = 3.0
    assert tm.dntp_conc == pytest.approx(3.0)


def test_base_pair_weights_tbl_non_table_equality() -> None:
    """Test BasePairWeightsTbl equality comparison with non-table objects."""
    from amplifyp.settings import DEFAULT_BASE_PAIR_WEIGHTS

    tbl = DEFAULT_BASE_PAIR_WEIGHTS
    assert tbl != "not_a_table"
    assert tbl != 42


def test_core_settings_dict_round_trip() -> None:
    """Test to_dict and from_dict round-tripping on core settings classes."""
    rep = ReplicationSettings(
        primability_cutoff=0.75,
        stability_cutoff=0.35,
        amplify4_compatibility_mode=True,
    )
    rep.base_pair_scores["A", "T"] = 42.0
    rep.match_weight[0] = 99.0
    rep.run_weights[1] = 111.0
    rep_round = ReplicationSettings.from_dict(rep.to_dict())
    assert rep_round == rep

    # String-typed GUI dictionary coercion
    rep_from_str = ReplicationSettings.from_dict(
        {
            "primability_cutoff": "0.65",
            "stability_cutoff": "bad",
            "amp4_compat": "true",
            "match_weight_default": "2",
            "match_weight_0": "45",
        }
    )
    assert rep_from_str.primability_cutoff == pytest.approx(0.65)
    assert rep_from_str.stability_cutoff == pytest.approx(0.4)
    assert rep_from_str.amplify4_compatibility_mode is True
    assert rep_from_str.match_weight[0] == pytest.approx(45.0)
    assert rep_from_str.match_weight[99] == pytest.approx(2.0)

    pd = PrimerDimerSettings(
        min_overlap=5, threshold=45.0, symbol_threshold=8.0
    )
    assert PrimerDimerSettings.from_dict(pd.to_dict()) == pd

    tm = TMSettings(
        dna_conc=100.0,
        dnap_conc=1.0,
        monovalent_salt_conc=75.0,
        divalent_salt_conc=2.0,
        dntp_conc=0.8,
    )
    tm.enthalpy["A", "A"] = 123.0
    tm.entropy["C", "C"] = 456.0
    assert TMSettings.from_dict(tm.to_dict()) == tm


def test_core_settings_partial_overrides_and_bad_coercions() -> None:
    """Partial overrides preserve defaults and reject bool/non-finite inputs."""
    from amplifyp.settings import (
        DEFAULT_BASE_PAIR_WEIGHTS,
        DEFAULT_MATCH_WEIGHTS,
        DEFAULT_PRIMABILITY_CUTOFF,
        DEFAULT_PRIMER_DIMER_OVERLAP,
    )

    rep_partial = ReplicationSettings.from_dict(
        {
            "bp_score_G_G": 90.0,
            "match_weight_0": 45.0,
            "primability_cutoff": True,
            "stability_cutoff": float("nan"),
        }
    )
    assert rep_partial.base_pair_scores["G", "G"] == pytest.approx(90.0)
    assert rep_partial.base_pair_scores["A", "A"] == pytest.approx(
        DEFAULT_BASE_PAIR_WEIGHTS["A", "A"]
    )
    assert rep_partial.match_weight[0] == pytest.approx(45.0)
    assert rep_partial.match_weight[1] == pytest.approx(
        DEFAULT_MATCH_WEIGHTS[1]
    )
    assert rep_partial.primability_cutoff == pytest.approx(
        DEFAULT_PRIMABILITY_CUTOFF
    )
    assert rep_partial.stability_cutoff == pytest.approx(0.4)

    rep_custom = ReplicationSettings()
    rep_custom.primer_dimer_scores["A", "T"] = 12.0
    assert ReplicationSettings.from_dict(rep_custom.to_dict()) == rep_custom

    pd_bool = PrimerDimerSettings.from_dict(
        {"pd_min_overlap": True, "pd_threshold": "inf"}
    )
    assert pd_bool.min_overlap == DEFAULT_PRIMER_DIMER_OVERLAP
    assert pd_bool.threshold == pytest.approx(60.0)

    rep_cleared = ReplicationSettings.from_dict({"match_weight_default": 2.0})
    assert rep_cleared.match_weight.default_weight == pytest.approx(2.0)
    assert rep_cleared.match_weight.overrides == {}

    # coerce_bool rejects NaN, and default matrix keys are omitted unless edited
    rep_nan = ReplicationSettings.from_dict({"amp4_compat": float("nan")})
    assert rep_nan.amplify4_compatibility_mode is False
    rep_bad_partial = ReplicationSettings.from_dict({"match_weight_0": "bad"})
    assert rep_bad_partial.match_weight[0] == pytest.approx(
        DEFAULT_MATCH_WEIGHTS[0]
    )
    default_rep_dict = ReplicationSettings().to_dict()
    assert not any(k.startswith("rep_pd_score_") for k in default_rep_dict)
    default_tm_dict = TMSettings().to_dict()
    assert not any(k.startswith("tm_enthalpy_") for k in default_tm_dict)
    assert not any(k.startswith("tm_entropy_") for k in default_tm_dict)

    # Non-string mapping keys are ignored safely
    assert (
        ReplicationSettings.from_dict({123: 456}) == ReplicationSettings()  # type: ignore[dict-item]
    )

    # GUISettings prunes omitted lengthwise weights and coerces
    # amplify4_compatibility_mode safely
    from amplifyp.gui.settings import GUISettings

    gui_settings = GUISettings()
    gui_settings.from_dict(
        {
            "match_weight_0": "45",
            "amplify4_compatibility_mode": float("nan"),
        }
    )
    assert gui_settings["match_weight_0"] == "45"
    assert "match_weight_1" not in gui_settings
    assert gui_settings["amp4_compat"] is False

    from amplifyp.settings import coerce_bool

    assert coerce_bool("", True) is True
    assert coerce_bool(1, False) is True
    assert coerce_bool(0, True) is False
    assert coerce_bool("maybe", True) is True

    pd_bad_int = PrimerDimerSettings.from_dict({"pd_min_overlap": "bad"})
    assert pd_bad_int.min_overlap == DEFAULT_PRIMER_DIMER_OVERLAP

    partial_default = BasePairWeightsTbl("AC", "AC", [[1.0, 2.0], [3.0, 4.0]])
    partial_default._col_map[ord("C")] = -1
    partial_default._col_map[ord("c")] = -1
    del partial_default._weight[("A", "C")]
    tbl_fallback = BasePairWeightsTbl.from_dict(
        {"tbl_A_A": 9.0}, "tbl", partial_default
    )
    assert tbl_fallback["A", "A"] == pytest.approx(9.0)
    assert tbl_fallback["A", "C"] == pytest.approx(0.0)


def test_overflow_coercion_and_gui_compat_normalisation() -> None:
    """Test overflow fallback and authoritative amp4_compat in GUISettings."""
    from amplifyp.gui.settings import GUISettings
    from amplifyp.settings import (
        DEFAULT_PRIMER_DIMER_OVERLAP,
        DEFAULT_PRIMER_DIMER_THRESHOLD,
    )

    pd_inf = PrimerDimerSettings.from_dict({"pd_min_overlap": float("inf")})
    assert pd_inf.min_overlap == DEFAULT_PRIMER_DIMER_OVERLAP

    pd_big = PrimerDimerSettings.from_dict({"pd_threshold": 10**400})
    assert pd_big.threshold == pytest.approx(DEFAULT_PRIMER_DIMER_THRESHOLD)

    gui_settings = GUISettings({"amplify4_compatibility_mode": True})
    assert "amplify4_compatibility_mode" not in gui_settings
    assert gui_settings["amp4_compat"] is True
    assert (
        gui_settings.get_replication_settings().amplify4_compatibility_mode
        is True
    )

    gui_settings["amp4_compat"] = False
    assert (
        gui_settings.get_replication_settings().amplify4_compatibility_mode
        is False
    )

    # Even if the legacy alias is set directly, amp4_compat is authoritative.
    gui_settings["amplify4_compatibility_mode"] = True
    assert (
        gui_settings.get_replication_settings().amplify4_compatibility_mode
        is False
    )

    # Float conversion in GUISettings.from_dict handles OverflowError safely.
    gui_settings.from_dict({"last_version_check_timestamp": 10**400})
    assert gui_settings["last_version_check_timestamp"] == pytest.approx(0.0)


def test_base_pair_weights_tbl_ascii_lut_invalidation_and_copy() -> None:
    """Verify ascii_lut caches correctly and invalidates on write or copy."""
    tbl = BasePairWeightsTbl("AT", "AT", [[1.0, -2.0], [-2.0, 1.0]])
    assert tbl.ascii_lut[ord("a")][ord("A")] == pytest.approx(1.0)
    assert tbl.ascii_lut[ord("Z")][ord("A")] == pytest.approx(0.0)

    tbl_copy = tbl.copy()
    tbl["A", "A"] = 3.5
    tbl["Z", "T"] = 4.25
    assert tbl.ascii_lut[ord("a")][ord("A")] == pytest.approx(3.5)
    assert tbl.ascii_lut[ord("Z")][ord("t")] == pytest.approx(4.25)
    assert tbl_copy.ascii_lut[ord("A")][ord("A")] == pytest.approx(1.0)
