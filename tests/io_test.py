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

"""Tests for amplifyp.io."""

from pathlib import Path

import pytest

from amplifyp.dna import DNAType
from amplifyp.io import (
    ProjectState,
    ReactionInput,
    get_duplicate_primer_indices,
    load_state_yaml,
    parse_primer_clipboard_text,
    parse_primers_from_text,
    parse_state_yaml,
    serialise_primers_to_tsv,
    serialise_state,
    validate_primer,
    validate_primers,
)
from amplifyp.settings import (
    PrimerDimerSettings,
    ReplicationSettings,
    TMSettings,
)


def test_validate_primer_and_case_sensitive_names() -> None:
    """Test single and batch validation with case-sensitive primer names."""
    assert validate_primer("", "", show_empty_errors=True) == (
        "Name cannot be empty",
        "Sequence cannot be empty",
    )
    assert validate_primer("P1", "INVALID_Z") == (
        None,
        "The DNA sequence contains invalid characters: I, L, Z, _",
    )
    assert validate_primer("P1", "ACGT") == (None, None)

    # Case-sensitive names: "p1" and "P1" are distinct; "P1" and "P1" collide.
    primers_distinct_case = [
        {"name": "p1", "seq": "ACGT", "active": True},
        {"name": "P1", "seq": "TGCA", "active": True},
    ]
    assert get_duplicate_primer_indices(primers_distinct_case) == set()
    assert validate_primers(primers_distinct_case) == [
        {"name": None, "seq": None},
        {"name": None, "seq": None},
    ]

    primers_dup = [
        {"name": "P1", "seq": "ACGT", "active": True},
        {"name": "P1", "seq": "acgt", "active": True},
        {"name": "P2", "seq": "TTTT", "active": False},
    ]
    assert get_duplicate_primer_indices(primers_dup) == {0, 1}
    errors = validate_primers(primers_dup)
    assert errors[0] == {
        "name": "Duplicate primer name",
        "seq": "Duplicate primer sequence",
    }
    assert errors[1] == {
        "name": "Duplicate primer name",
        "seq": "Duplicate primer sequence",
    }
    assert errors[2] == {"name": None, "seq": None}


def test_primer_tabular_and_clipboard_io() -> None:
    """Test CSV/TSV primer parsing, TSV serialisation, and clipboard parsing."""
    raw = (
        "# comment\n"
        "\n"
        "ACGT\tFwd1\n"
        "TGCA,Rev1,Extra Note\n"
        "INVALIDZ,BadPrimer\n"
        "singlecol\n"
    )
    parsed = parse_primers_from_text(raw)
    assert parsed == [
        {"name": "Fwd1", "seq": "ACGT", "active": True},
        {"name": "Rev1 - Extra Note", "seq": "TGCA", "active": True},
        {"name": "BadPrimer", "seq": "INVALIDZ", "active": False},
    ]

    tsv = serialise_primers_to_tsv(parsed)
    assert "ACGT\tFwd1" in tsv
    assert "TGCA\tRev1 - Extra Note" in tsv

    clip = parse_primer_clipboard_text("P1\tACGT\nP2,TGCA\nRYMKSW\nJustAName\n")
    assert clip == [
        {"name": "P1", "seq": "ACGT", "active": False},
        {"name": "P2", "seq": "TGCA", "active": False},
        {"name": "", "seq": "RYMKSW", "active": False},
        {"name": "JustAName", "seq": "", "active": False},
    ]


def test_reaction_input_and_pcr_conversion() -> None:
    """Test ReactionInput state helpers and core PCR/dimer builders."""
    tpl = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    ri = ReactionInput(
        template=tpl,
        template_circular=True,
        primers=[
            {"name": "F", "seq": "ACGTACGT", "active": True},
            {"name": "R", "seq": "ACGTACGT", "active": False},
        ],
    )
    assert len(ri.get_active_primers()) == 1
    assert ri.to_template_dna().type == DNAType.CIRCULAR

    ri.activate_only_primers([("F_new", "ACGTACGT"), ("R_new", "TGCATGCA")])
    active = ri.get_active_primers()
    assert [(p["name"], p["seq"]) for p in active] == [
        ("F_new", "ACGTACGT"),
        ("R_new", "TGCATGCA"),
    ]

    pcr = ri.to_pcr()
    assert len(pcr) == 2
    gen = ri.to_dimer_generator()
    assert len(gen.primers) == 2


def test_project_state_yaml_round_trip(tmp_path: Path) -> None:
    """Test serialise_state, parse_state_yaml, and load_state_yaml."""
    simple_yaml = Path("tests/examples/save_states/simple.yaml")
    state = load_state_yaml(simple_yaml)
    assert len(state.input.primers) == 4
    assert state.input.template_circular is False

    dumped = state.to_yaml(include_settings=True)
    reloaded = parse_state_yaml(dumped)
    assert reloaded.input.to_dict() == state.input.to_dict()
    assert reloaded.replication_settings == ReplicationSettings()
    assert reloaded.primer_dimer_settings == PrimerDimerSettings()
    assert reloaded.tm_settings == TMSettings()

    yaml_text = serialise_state({"multiline": "line1\nline2", "plain": "abc"})
    assert "multiline: |-" in yaml_text or "multiline: |" in yaml_text

    bad_path = tmp_path / "bad.yaml"
    bad_path.write_text("- not a mapping\n", encoding="utf-8")
    with pytest.raises(ValueError, match="expected a mapping"):
        load_state_yaml(bad_path)

    # Default ProjectState construction and null-safe ReactionInput.from_dict
    default_state = ProjectState(input=ReactionInput())
    assert "input" in default_state.to_dict()

    ri_null = ReactionInput(
        template="ACGT", primers=[{"name": "P", "seq": "ACGT"}]
    )
    ri_null.from_dict({"template": None, "primers": None})
    assert ri_null.template == ""
    assert ri_null.primers == []

    # Case-insensitive matching in activate_only_primers
    ri_case = ReactionInput(
        primers=[{"name": "Old", "seq": "acgtacgt", "active": False}]
    )
    ri_case.activate_only_primers([("Renamed", "ACGTACGT")])
    assert len(ri_case.primers) == 1
    assert ri_case.primers[0]["name"] == "Renamed"
    assert ri_case.primers[0]["active"] is True

    # Whitespace stripped from CSV/TSV sequence column
    spaced = parse_primers_from_text("ACGT TGCA,Spaced\n")
    assert spaced == [{"name": "Spaced", "seq": "ACGTTGCA", "active": True}]

    # Coerces string "false" for template_circular and validates bad inputs
    ri_bool = ReactionInput()
    ri_bool.from_dict({"template_circular": "false"})
    assert ri_bool.template_circular is False
    with pytest.raises(ValueError, match="sequence"):
        ri_bool.from_dict({"primers": "ACGT"})
    with pytest.raises(ValueError, match="mapping"):
        ri_bool.from_dict({"primers": [None]})
    with pytest.raises(ValueError, match="'input' must be a mapping"):
        ProjectState.from_dict({"input": []})
    with pytest.raises(ValueError, match="'settings' must be a mapping"):
        ProjectState.from_dict({"settings": []})

    # ProjectState keeps amp4_compat synced with amplify4_compatibility_mode and
    # GUISettings honours amplify4_compatibility_mode
    from amplifyp.gui.settings import GUISettings

    synced_state = ProjectState.from_dict({"settings": {"amp4_compat": False}})
    synced_state.replication_settings.amplify4_compatibility_mode = True
    gui_settings = GUISettings()
    gui_settings.from_dict(synced_state.to_dict()["settings"])
    assert gui_settings.get_replication_settings().amplify4_compatibility_mode

    core_only_state = ProjectState(
        replication_settings=ReplicationSettings(
            amplify4_compatibility_mode=True
        )
    )
    gui_settings_core = GUISettings()
    gui_settings_core.from_dict(core_only_state.to_dict()["settings"])
    assert (
        gui_settings_core.get_replication_settings().amplify4_compatibility_mode
        is True
    )

    # StrEnum values serialise safely and round-trip through parse_state_yaml
    from amplifyp.melting import TmMethod

    enum_yaml = serialise_state(
        {"settings": {"tm_method": TmMethod.LANDER_AMPLIFY4}}
    )
    assert "!!python/" not in enum_yaml
    reloaded_enum = parse_state_yaml(enum_yaml)
    assert reloaded_enum.raw_settings["tm_method"] == "Lander / Amplify 4"

    # List-valued primer seq, empty candidate in activate_only_primers,
    # non-mapping ReactionInput.from_dict, and non-string YAML settings keys
    ri_list = ReactionInput()
    ri_list.from_dict(
        {"primers": [{"name": "Split", "seq": ["AC", "GT"], "active": True}]}
    )
    assert ri_list.primers == [{"name": "Split", "seq": "ACGT", "active": True}]
    ri_list.activate_only_primers([("Empty", "   "), ("Split", "ACGT")])
    assert len(ri_list.get_active_primers()) == 1
    with pytest.raises(ValueError, match="expected a mapping"):
        ri_list.from_dict([])  # type: ignore[arg-type]

    non_str_state = ProjectState.from_dict({"settings": {123: 456}})
    assert non_str_state.replication_settings == ReplicationSettings()

    # Escaped whitespace counts as empty in validate_primer
    assert validate_primer("P1", "\\n\\t", show_empty_errors=True) == (
        None,
        "Sequence cannot be empty",
    )

    # serialise_primers_to_tsv cleans multiline sequences and uses LF endings
    long_seq = "ACGT" * 25
    tsv_wrapped = serialise_primers_to_tsv(
        [{"name": "Long\nPrimer", "seq": f"{long_seq[:80]}\n{long_seq[80:]}"}]
    )
    assert "\r" not in tsv_wrapped
    assert parse_primers_from_text(tsv_wrapped) == [
        {"name": "Long Primer", "seq": long_seq, "active": True}
    ]
