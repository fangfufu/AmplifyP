# Copyright (C) 2026 Fufu Fang
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

"""Tests for GUI state saving and loading."""

from unittest.mock import MagicMock

import flet as ft
import yaml

from amplifyp.gui.views.input_view import InputView
from amplifyp.gui.views.settings_view import SettingsView


def test_gui_state_save_load() -> None:
    """Test saving and loading GUI state."""
    mock_page = MagicMock(spec=ft.Page)

    # 1. Setup InputView with data
    input_view = InputView(mock_page)

    # Set template sequence (long enough to trigger multiline formatting)
    template_seq = "ATGC" * 30
    input_view.template_sequence.value = template_seq
    input_view.template_circular.value = True

    # Add primers
    input_view.primer_name_input.value = "Primer1"
    input_view.primer_seq_input.value = "ATGC"
    input_view.add_primer_clicked(MagicMock())

    input_view.primer_name_input.value = "Primer2"
    input_view.primer_seq_input.value = "CGTA"
    # Simulate unchecked primer
    # We can't easily simulate UI interaction to uncheck, so we modify
    # the control directly.
    # add_primer_clicked adds a checked primer by default.
    input_view.add_primer_clicked(MagicMock())

    # Find the second primer checkbox and uncheck it
    # accessing the last added primer
    last_tile = input_view.primers_list.controls[-1]
    if isinstance(last_tile, ft.ListTile) and isinstance(
        last_tile.leading, ft.Checkbox
    ):
        last_tile.leading.value = False

    # 2. Setup SettingsView with modified settings
    settings_view = SettingsView(mock_page)
    settings_view.set_primability_cutoff.value = "0.9"
    settings_view.set_amp4_compat.value = True
    settings_view.set_tm_dna_conc.value = "100.0"

    # 3. Capture State
    input_state = input_view.get_state()
    settings_state = settings_view.get_state()

    # Combine state as main.py does
    full_state = input_state.copy()
    full_state["settings"] = settings_state

    # 4. Serialize (replicating main.py logic)
    def multiline_presenter(dumper: yaml.Dumper, data: str) -> yaml.ScalarNode:
        if "\n" in data:
            return dumper.represent_scalar(
                "tag:yaml.org,2002:str", data, style="|"
            )
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    yaml.add_representer(str, multiline_presenter)
    yaml_str = yaml.dump(full_state, sort_keys=False)

    # 5. Deserialize
    loaded_state = yaml.safe_load(yaml_str)

    # 6. Restore to new views
    new_input_view = InputView(mock_page)
    new_settings_view = SettingsView(mock_page)

    new_input_view.set_state(loaded_state)
    new_settings_view.set_state(loaded_state)

    # 7. Assertions

    # Verify InputView
    # Template sequence should match (ignoring newlines inserted by formatting)
    assert (
        new_input_view.template_sequence.value.replace("\n", "") == template_seq
    )
    assert new_input_view.template_circular.value

    # Verify Primers
    # get_primers() only returns active primers, so we should check controls
    # directly or get_all_primers_state.
    # But let's check what get_primers returns first
    active_primers = new_input_view.get_primers()
    assert len(active_primers) == 1
    assert active_primers[0]["name"] == "Primer1"
    assert active_primers[0]["seq"] == "ATGC"

    # Check all primers including inactive ones
    all_primers = new_input_view.get_all_primers_state()
    assert len(all_primers) == 2

    p1 = next(p for p in all_primers if p["name"] == "Primer1")
    assert p1["seq"].replace("\n", "") == "ATGC"
    assert p1["active"]

    p2 = next(p for p in all_primers if p["name"] == "Primer2")
    assert p2["seq"].replace("\n", "") == "CGTA"
    assert not p2["active"]

    # Verify SettingsView
    assert new_settings_view.set_primability_cutoff.value == "0.9"
    assert new_settings_view.set_amp4_compat.value
    assert new_settings_view.set_tm_dna_conc.value == "100.0"
    # Check a default value wasn't changed
    assert new_settings_view.set_stability_cutoff.value == "0.4"
