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

from amplifyp.gui.views.input import InputView
from amplifyp.gui.views.settings_view import SettingsView


def test_gui_state_save_load() -> None:
    """Test saving and loading GUI state."""
    mock_page = MagicMock(spec=ft.Page)

    # 1. Setup InputView with data
    input_view = InputView(mock_page)

    # Set template sequence (long enough to trigger multiline formatting)
    template_seq = "ATGC" * 30
    input_view.input_data.template = template_seq
    input_view.input_data.template_circular = True

    # Add primers
    input_view.input_data.primers = [
        {"name": "Primer1", "seq": "ATGC", "active": True},
        {"name": "Primer2", "seq": "CGTA", "active": False},
    ]
    input_view.update_ui()

    # 2. Setup SettingsView with modified settings
    settings_view = SettingsView(mock_page)
    settings_view.set_primability_cutoff.value = "0.9"
    settings_view.set_amp4_compat.value = True
    settings_view.set_improved_visualisation.value = True
    settings_view.set_tm_dna_conc.value = "100.0"
    settings_view.set_tm_method.value = "Lander / Amplify 4"
    settings_view.set_font_family.value = "Courier New"
    settings_view.appearance_tile.set_color_scheme.value = (
        "Dark (Colour Deficient Friendly)"
    )
    settings_view.settings_map["bp_score_G_G"].value = "99.0"
    settings_view.settings_map["pd_score_G_G"].value = "99.0"

    # 3. Capture State - now as separate input + settings dicts
    input_state = input_view.get_state()
    settings_state = settings_view.get_state()

    # Combine state as app.py does: {"input": ..., "settings": ...}
    full_state = {
        "input": input_state,
        "settings": settings_state,
    }

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

    new_input_view.set_state(loaded_state["input"])
    new_settings_view.set_state(loaded_state["settings"])

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
    assert new_settings_view.set_improved_visualisation.value
    assert new_settings_view.set_tm_dna_conc.value == "100.0"
    assert new_settings_view.set_tm_method.value == "Lander / Amplify 4"
    assert new_settings_view.set_font_family.value == "Courier New"
    assert (
        new_settings_view.appearance_tile.set_color_scheme.value
        == "Dark (Colour Deficient Friendly)"
    )
    assert new_settings_view.settings["dark_mode"] is True
    assert new_settings_view.settings["color_deficient"] is True
    assert new_settings_view.settings["improved_visualisation"] is True
    assert new_settings_view.settings_map["bp_score_G_G"].value == "99.0"
    assert new_settings_view.settings_map["pd_score_G_G"].value == "99.0"
    # Check a default value wasn't changed
    assert new_settings_view.set_stability_cutoff.value == "0.4"

    # Reset GUIColors to avoid test contamination
    from amplifyp.gui.colors import GUIColors

    GUIColors.color_deficient_mode = False
    GUIColors.dark_mode = False


def test_settings_view_buttons() -> None:
    """Test Apply and Reset to Default buttons in SettingsView."""
    mock_page = MagicMock(spec=ft.Page)

    apply_called = False
    reset_called = False

    def on_apply_callback(e: ft.ControlEvent) -> None:
        nonlocal apply_called
        apply_called = True

    def on_reset_callback(e: ft.ControlEvent) -> None:
        nonlocal reset_called
        reset_called = True

    settings_view = SettingsView(
        mock_page,
        on_apply=on_apply_callback,
        on_reset=on_reset_callback,
    )

    # Change some values
    settings_view.set_primability_cutoff.value = "0.95"
    settings_view.set_amp4_compat.value = True
    settings_view.set_improved_visualisation.value = True
    settings_view.settings_map["bp_score_G_G"].value = "50.0"
    settings_view.settings_map["pd_score_G_G"].value = "50.0"

    # Find the Row containing the Apply and Reset buttons
    buttons_row = settings_view.controls[-1]
    assert isinstance(buttons_row, ft.Row)
    apply_btn = buttons_row.controls[0]
    reset_btn = buttons_row.controls[1]

    assert apply_btn.content == "Apply"
    assert reset_btn.content == "Reset to Default"

    # Trigger Apply
    apply_btn.on_click(MagicMock(spec=ft.ControlEvent))
    assert apply_called
    assert settings_view.settings["primability_cutoff"] == "0.95"
    assert settings_view.settings["amp4_compat"] is True
    assert settings_view.settings["improved_visualisation"] is True
    assert settings_view.settings["bp_score_G_G"] == "50.0"
    assert settings_view.settings["pd_score_G_G"] == "50.0"

    # Trigger Reset
    reset_btn.on_click(MagicMock(spec=ft.ControlEvent))
    assert reset_called
    assert (
        settings_view.settings["tm_method"]
        == "SantaLucia 1998 / Owczarzy 2008 (Default)"
    )
    assert settings_view.settings["bp_score_G_G"] == "100"
    assert settings_view.settings["pd_score_G_G"] == "-20"
    assert settings_view.settings["improved_visualisation"] is True
    # Controls should be updated too
    assert settings_view.set_primability_cutoff.value == "0.8"
    assert settings_view.set_amp4_compat.value is False
    assert settings_view.set_improved_visualisation.value is True
    assert (
        settings_view.set_tm_method.value
        == "SantaLucia 1998 / Owczarzy 2008 (Default)"
    )
    assert settings_view.settings_map["bp_score_G_G"].value == "100"
    assert settings_view.settings_map["pd_score_G_G"].value == "-20"


def test_color_deficient_mode_switching() -> None:
    """Test toggling color deficient setting shifts GUIColors."""
    from amplifyp.gui.colors import GUIColors
    from amplifyp.gui.settings import GUISettings

    settings = GUISettings()
    # 1. Initially false/standard
    assert settings["color_deficient"] is False
    assert GUIColors.color_deficient_mode is False
    standard_success = GUIColors.SUCCESS_GREEN
    standard_error = GUIColors.ERROR_RED
    standard_fwd = GUIColors.FWD_PRIMER
    standard_rev = GUIColors.REV_PRIMER

    # 2. Toggle setting to True
    settings["color_deficient"] = True
    assert GUIColors.color_deficient_mode is True
    assert GUIColors.SUCCESS_GREEN != standard_success
    assert GUIColors.ERROR_RED != standard_error
    assert GUIColors.FWD_PRIMER != standard_fwd
    assert GUIColors.REV_PRIMER != standard_rev

    # 3. Toggle back
    settings["color_deficient"] = False
    assert GUIColors.color_deficient_mode is False
    assert GUIColors.SUCCESS_GREEN == standard_success
    assert GUIColors.ERROR_RED == standard_error
    assert GUIColors.FWD_PRIMER == standard_fwd
    assert GUIColors.REV_PRIMER == standard_rev


def test_dark_mode_switching() -> None:
    """Test toggling dark mode setting shifts GUIColors."""
    from amplifyp.gui.colors import GUIColors
    from amplifyp.gui.settings import GUISettings

    settings = GUISettings()
    # 1. Initially false/standard
    assert settings["dark_mode"] == "system"
    assert GUIColors.dark_mode is False
    standard_selected_bg = GUIColors.SELECTED_ROW_BG
    standard_info_header = GUIColors.INFO_HEADER_BG
    standard_fwd = GUIColors.FWD_PRIMER
    standard_rev = GUIColors.REV_PRIMER

    # 2. Toggle setting to True
    settings["dark_mode"] = True
    assert GUIColors.dark_mode is True
    assert GUIColors.SELECTED_ROW_BG != standard_selected_bg
    assert GUIColors.INFO_HEADER_BG != standard_info_header
    assert GUIColors.FWD_PRIMER != standard_fwd
    assert GUIColors.REV_PRIMER != standard_rev

    # 3. Toggle back
    settings["dark_mode"] = False
    assert GUIColors.dark_mode is False
    assert GUIColors.SELECTED_ROW_BG == standard_selected_bg
    assert GUIColors.INFO_HEADER_BG == standard_info_header
    assert GUIColors.FWD_PRIMER == standard_fwd
    assert GUIColors.REV_PRIMER == standard_rev


def test_system_theme_saving_loading() -> None:
    """Test saving and loading system theme settings."""
    mock_page = MagicMock(spec=ft.Page)
    settings_view = SettingsView(mock_page)

    # 1. Select System dropdown
    settings_view.appearance_tile.set_color_scheme.value = "System"
    settings_view.sync_to_state()

    assert settings_view.settings["dark_mode"] == "system"
    assert settings_view.settings["color_deficient"] is False

    # 2. Serialize / deserialize settings state
    settings_state = settings_view.get_state()
    new_settings_view = SettingsView(mock_page)
    new_settings_view.set_state(settings_state)

    # 3. Assertions
    assert new_settings_view.appearance_tile.set_color_scheme.value == "System"
    assert new_settings_view.settings["dark_mode"] == "system"
    assert new_settings_view.settings["color_deficient"] is False

    # 4. Now with Colour Deficient Friendly System option
    settings_view.appearance_tile.set_color_scheme.value = (
        "System (Colour Deficient Friendly)"
    )
    settings_view.sync_to_state()

    assert settings_view.settings["dark_mode"] == "system"
    assert settings_view.settings["color_deficient"] is True

    settings_state_2 = settings_view.get_state()
    new_settings_view_2 = SettingsView(mock_page)
    new_settings_view_2.set_state(settings_state_2)

    assert (
        new_settings_view_2.appearance_tile.set_color_scheme.value
        == "System (Colour Deficient Friendly)"
    )
    assert new_settings_view_2.settings["dark_mode"] == "system"
    assert new_settings_view_2.settings["color_deficient"] is True

    # Reset GUIColors to avoid test contamination
    from amplifyp.gui.colors import GUIColors

    GUIColors.color_deficient_mode = False
    GUIColors.dark_mode = False


def test_simple_state_font_sizes_are_integers() -> None:
    """Test loading simple.yaml parses font sizes as integers."""
    import os

    from amplifyp.gui.settings import GUISettings

    yaml_path = os.path.join(
        os.path.dirname(__file__), "examples", "save_states", "simple.yaml"
    )
    with open(yaml_path, encoding="utf-8") as f:
        parsed_state = yaml.safe_load(f)

    settings = GUISettings()
    settings.from_dict(parsed_state["settings"])

    # Assert that all font size settings are actual integers
    font_size_keys = [
        "font_size_map_baseline",
        "font_size_map_primer",
        "font_size_map_amplicon",
        "font_size_header",
        "font_size_subheader",
        "font_size_body",
        "font_size_small",
        "font_size_default",
        "font_size_micro",
        "font_size_table_header",
    ]
    for key in font_size_keys:
        val = settings.get(key)
        assert isinstance(val, int), (
            f"{key} should be an int, but got {type(val)}"
        )

    # Test manual assignment using __setitem__ with string values
    settings["font_size_map_baseline"] = "20"
    assert isinstance(settings.get("font_size_map_baseline"), int)
    assert settings.get("font_size_map_baseline") == 20
