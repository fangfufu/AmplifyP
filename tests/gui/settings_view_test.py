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

"""Tests for Settings View and sub-tiles."""

import asyncio
import os
from typing import Any
from unittest.mock import AsyncMock, MagicMock, patch

import flet as ft
import yaml

from amplifyp.gui.settings import GUISettings
from amplifyp.gui.views.settings import SettingsView
from amplifyp.gui.views.settings.diagnostics_tile import DiagnosticsTile
from amplifyp.gui.views.settings.general_tile import GeneralTile
from amplifyp.gui.views.settings.input_tile import InputTile


def test_settings_view_properties_and_methods() -> None:
    """Test all properties, sync, reset, and state methods in SettingsView."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    mock_on_change = MagicMock()
    mock_on_reset = MagicMock()

    settings = GUISettings()
    view = SettingsView(
        page=mock_page,
        settings=settings,
        on_change=mock_on_change,
        on_reset=mock_on_reset,
    )

    # 1. Properties and tile order
    assert hasattr(view, "general_tile")
    assert hasattr(view, "input_tile")
    assert view.controls[0] == view.general_tile
    assert isinstance(view.controls[1], ft.Divider)
    assert view.controls[2] == view.input_tile

    assert view.general_tile.title.value == "General"
    assert view.input_tile.title.value == "Input"
    assert view.replication_tile.title.value == "PCR"
    assert view.dimer_tile.title.value == "Primer Dimers"
    assert view.designer_2d_tile.title.value == "Designer 2D"
    assert view.diagnostics_tile.title.value == "Diagnostics"

    assert view.set_primability_cutoff is not None
    assert view.set_stability_cutoff is not None
    assert view.set_amp4_compat is not None
    assert view.set_tm_dna_conc is not None
    assert view.set_tm_dnap_conc is not None
    assert view.set_tm_mono_salt is not None
    assert view.set_tm_div_salt is not None
    assert view.set_tm_dNTP_conc is not None
    assert view.set_pd_min_overlap is not None
    assert view.set_pd_threshold is not None
    assert view.set_font_family is not None
    assert view.set_colour_scheme is not None
    assert view.set_sequence_separator is not None
    assert view.set_sequence_separator.value == "Dash ('-')"
    assert [opt.key for opt in view.set_sequence_separator.options] == [
        "Space (' ')",
        "Dash ('-')",
    ]
    assert view.set_colour_deficient is not None
    assert view.set_improved_visualisation is not None
    assert view.set_show_primer_temperature is not None
    assert view.set_tm_colour_scheme is not None
    assert view.set_designer_2d_colour_scheme is not None
    assert view.set_designer_2d_show_rev_fwd is not None
    assert view.set_version_checking_frequency is not None
    assert view.set_version_checking_frequency.value == "Once per Day"
    assert settings["version_checking_frequency"] == "Once per Day"
    assert view.set_auto_reload_on_startup is not None
    assert view.set_primer_info_panel_position is not None
    assert view.set_primer_info_panel_fixed_height is not None
    assert view.set_auto_activate_new_valid_primer is not None

    # 2. sync_to_state with invalid log_rotation_max_bytes
    view.diagnostics_tile.log_rotation_max_bytes.value = "invalid_mb"
    view.sync_to_state()
    assert settings["log_rotation_max_bytes"] == 5242880

    # 3. _on_change_handler with control match and callback invocation
    mock_event = MagicMock()
    mock_event.control = view.set_tm_method
    view.set_tm_method.value = "SantaLucia 1998 / Owczarzy 2008 (Default)"
    view._on_change_handler(mock_event)
    assert mock_on_change.called

    # 4. _on_reset_handler
    settings["designer_1d_show_tm"] = False
    settings["designer_1d_show_pct_at"] = True
    settings["pcr_amplicon_ranking"] = "Quality score"
    settings["pcr_relative_q_score"] = False
    settings["pcr_primer_label_orientation"] = "Vertical"
    view._on_reset_handler(mock_event)
    assert mock_on_reset.called
    assert settings["designer_1d_show_tm"] is True
    assert settings["designer_1d_show_pct_at"] is False
    assert settings["pcr_amplicon_ranking"] == "Position, then length"
    assert settings["pcr_relative_q_score"] is True
    assert settings["pcr_primer_label_orientation"] == "Horizontal"
    assert view.set_designer_1d_show_tm.value is True
    assert view.set_designer_1d_show_pct_at.value is False
    assert view.set_pcr_amplicon_ranking.value == "Position, then length"
    assert view.set_pcr_relative_q_score.value is True
    assert view.set_pcr_primer_label_orientation.value == "Horizontal"

    # 5. get_replication_settings, get_state, set_state
    repl_settings = view.get_replication_settings()
    assert repl_settings is not None

    state = view.get_state()
    assert isinstance(state, dict)
    state["font_family"] = "Courier New"
    view.set_state(state)
    assert settings["font_family"] == "Courier New"


def test_diagnostics_tile_all_branches(tmp_path: Any) -> None:
    """Test all branches of DiagnosticsTile."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    mock_page.services = []
    mock_on_change = MagicMock()

    # Initialise with non-standard max_bytes (15MB)
    settings = GUISettings()
    settings["log_rotation_max_bytes"] = 15 * 1024 * 1024
    log_file = tmp_path / "custom_app.log"
    settings["log_file_path"] = str(log_file)

    tile = DiagnosticsTile(
        page=mock_page,
        settings=settings,
        settings_map={},
        on_change_handler=mock_on_change,
        header_size=18,
    )
    assert "15" in [opt.key for opt in tile.log_rotation_max_bytes.options]
    assert any(opt.key == str(log_file) for opt in tile.log_file_path.options)

    # 1. _on_file_logging_change
    tile._on_file_logging_change(MagicMock())
    assert mock_on_change.called

    # 2. _on_log_file_path_change with regular path
    tile.log_file_path.value = "(Default)"
    asyncio.run(tile._on_log_file_path_change(MagicMock()))

    # 3. _on_log_file_path_change with "Select folder" - successful selection
    tile.log_file_path.value = "Select folder"
    new_dir = str(tmp_path / "new_log_dir")
    expected_log_path = os.path.join(new_dir, "app.log")
    with patch.object(
        ft.FilePicker,
        "get_directory_path",
        new=AsyncMock(return_value=new_dir),
    ):
        asyncio.run(tile._on_log_file_path_change(MagicMock()))
        assert settings["log_file_path"] == expected_log_path

    # 4. _on_log_file_path_change with "Select folder" - cancelled selection
    tile.log_file_path.value = "Select folder"
    with patch.object(
        ft.FilePicker, "get_directory_path", new=AsyncMock(return_value=None)
    ):
        asyncio.run(tile._on_log_file_path_change(MagicMock()))
        assert tile.log_file_path.value == expected_log_path

    # 5. update_ui with invalid max_bytes and missing mb option
    settings["log_rotation_max_bytes"] = "invalid_bytes"
    tile.update_ui()
    assert tile.log_rotation_max_bytes.value == "5"

    settings["log_rotation_max_bytes"] = 30 * 1024 * 1024
    tile.update_ui()
    assert tile.log_rotation_max_bytes.value == "30"

    # 6. Web mode
    mock_web_page = MagicMock(spec=ft.Page)
    mock_web_page.web = True
    web_tile = DiagnosticsTile(
        page=mock_web_page,
        settings=settings,
        settings_map={},
        on_change_handler=mock_on_change,
        header_size=18,
    )
    assert web_tile.log_file_path.visible is False
    assert web_tile.log_rotation_enabled.visible is False


def test_general_tile_all_branches() -> None:
    """Test all branches of GeneralTile including update checker and IO."""
    mock_page = MagicMock(spec=ft.Page)
    mock_on_change = MagicMock()
    mock_sync = MagicMock()
    mock_update_ui = MagicMock()
    mock_update_found = MagicMock()

    settings = GUISettings()
    tile = GeneralTile(
        page=mock_page,
        settings=settings,
        settings_map={},
        on_change_handler=mock_on_change,
        header_size=18,
        font_size_default=14,
        sync_to_state_callback=mock_sync,
        update_ui_callback=mock_update_ui,
        on_update_found=mock_update_found,
    )

    # 1. Properties
    assert tile.set_colour_deficient is not None
    assert tile.set_auto_reload_on_startup is not None
    assert tile.set_sequence_separator is not None
    assert tile.set_sequence_separator.value == "Dash ('-')"

    # 2. _on_manual_check_click
    captured_task = None

    def mock_run_task(task: Any) -> None:
        nonlocal captured_task
        captured_task = task

    mock_page.run_task = mock_run_task
    tile._on_manual_check_click(MagicMock())
    assert captured_task == tile.perform_manual_check

    # 3. perform_manual_check: network error (None)
    with patch(
        "amplifyp.gui.utils.system.fetch_latest_release_version",
        return_value=None,
    ):
        asyncio.run(tile.perform_manual_check())
        assert "Could not check for updates" in tile.status_text.value

    # 4. perform_manual_check: newer version found
    with (
        patch(
            "amplifyp.gui.utils.system.fetch_latest_release_version",
            return_value="v99.0.0",
        ),
        patch("amplifyp.gui.utils.system.is_newer_version", return_value=True),
    ):
        asyncio.run(tile.perform_manual_check())
        assert "New version v99.0.0 is available!" in tile.status_text.value
        mock_update_found.assert_called_with("v99.0.0")

    # 5. perform_manual_check: up to date
    with (
        patch(
            "amplifyp.gui.utils.system.fetch_latest_release_version",
            return_value="v0.1.0",
        ),
        patch("amplifyp.gui.utils.system.is_newer_version", return_value=False),
    ):
        asyncio.run(tile.perform_manual_check())
        assert "AmplifyP is up to date." in tile.status_text.value

    # 6. _save_settings_async: filepicker open, success, and error handling
    tile.filepicker_open = True
    asyncio.run(tile._save_settings_async(MagicMock()))
    tile.filepicker_open = False

    with patch(
        "amplifyp.gui.views.settings.general_tile.save_and_write_file",
        new=AsyncMock(return_value=True),
    ):
        asyncio.run(tile._save_settings_async(MagicMock()))

    with (
        patch(
            "amplifyp.gui.views.settings.general_tile.save_and_write_file",
            new=AsyncMock(side_effect=OSError("Save failed")),
        ),
        patch.object(tile.notification_helper, "show_message") as mock_msg,
    ):
        asyncio.run(tile._save_settings_async(MagicMock()))
        assert any(
            "Error saving settings" in str(arg) for arg in mock_msg.call_args[0]
        )

    # 7. _load_settings_async: filepicker open, cancelled, invalid, and errors
    tile.filepicker_open = True
    asyncio.run(tile._load_settings_async(MagicMock()))
    tile.filepicker_open = False

    # Cancelled
    with patch(
        "amplifyp.gui.views.settings.general_tile.pick_and_read_file",
        new=AsyncMock(return_value=None),
    ):
        asyncio.run(tile._load_settings_async(MagicMock()))

    # Not a dict
    with (
        patch(
            "amplifyp.gui.views.settings.general_tile.pick_and_read_file",
            new=AsyncMock(return_value="[1, 2, 3]"),
        ),
        patch.object(tile.notification_helper, "show_message") as mock_msg,
    ):
        asyncio.run(tile._load_settings_async(MagicMock()))
        assert any(
            "Invalid settings file format" in str(arg)
            for arg in mock_msg.call_args[0]
        )

    # settings key not a dict
    with (
        patch(
            "amplifyp.gui.views.settings.general_tile.pick_and_read_file",
            new=AsyncMock(return_value="settings: 'invalid'"),
        ),
        patch.object(tile.notification_helper, "show_message") as mock_msg,
    ):
        asyncio.run(tile._load_settings_async(MagicMock()))
        assert any(
            "Invalid settings data format" in str(arg)
            for arg in mock_msg.call_args[0]
        )

    # Exception during load
    with (
        patch(
            "amplifyp.gui.views.settings.general_tile.pick_and_read_file",
            new=AsyncMock(return_value="valid: yaml"),
        ),
        patch("yaml.safe_load", side_effect=yaml.YAMLError("YAML fail")),
        patch.object(tile.notification_helper, "show_message") as mock_msg,
    ):
        asyncio.run(tile._load_settings_async(MagicMock()))
        assert any(
            "Error loading settings" in str(arg)
            for arg in mock_msg.call_args[0]
        )

    # Successful load with nested settings dict
    with (
        patch(
            "amplifyp.gui.views.settings.general_tile.pick_and_read_file",
            new=AsyncMock(
                return_value="settings:\n  font_family: 'Courier New'"
            ),
        ),
        patch.object(tile.notification_helper, "show_message") as mock_msg,
    ):
        asyncio.run(tile._load_settings_async(MagicMock()))
        assert settings["font_family"] == "Courier New"
        assert any(
            "Settings loaded successfully" in str(arg)
            for arg in mock_msg.call_args[0]
        )

    # Successful load with flat dict (no settings key)
    with (
        patch(
            "amplifyp.gui.views.settings.general_tile.pick_and_read_file",
            new=AsyncMock(return_value="font_family: 'Roboto Mono'"),
        ),
        patch.object(tile.notification_helper, "show_message") as mock_msg,
    ):
        asyncio.run(tile._load_settings_async(MagicMock()))
        assert settings["font_family"] == "Roboto Mono"

    # 8. Colour scheme dropdown and syncing (all 6 branches)
    schemes = [
        ("Dark", True, False),
        ("Dark (Colour Deficient Friendly)", True, True),
        ("Light (Colour Deficient Friendly)", False, True),
        ("System", "system", False),
        ("System (Colour Deficient Friendly)", "system", True),
        ("Light", False, False),
    ]
    for scheme_val, expected_dark, expected_def in schemes:
        tile.set_colour_scheme.value = scheme_val
        tile._on_colour_scheme_change(MagicMock())
        assert settings["dark_mode"] == expected_dark
        assert settings["colour_deficient"] == expected_def
        tile.update_colour_scheme_dropdown()
        assert tile.set_colour_scheme.value == scheme_val

    # 9. GeneralTile appearance controls and update_ui
    assert tile.set_font_family is not None
    assert tile.set_colour_scheme is not None
    assert tile.set_colour_deficient is not None
    assert tile.set_sequence_separator is not None
    settings["sequence_separator"] = "Space (' ')"
    settings["dark_mode"] = True
    settings["colour_deficient"] = False
    tile.update_ui()
    assert tile.set_sequence_separator.value == "Space (' ')"
    assert tile.set_colour_scheme.value == "Dark"

    from amplifyp.gui.colours import GUIColours

    tile._status_state = "error"
    tile.update_ui()
    assert tile.status_text.color == GUIColours.ERROR_RED

    tile._status_state = "checking"
    tile.update_ui()
    assert tile.status_text.color == GUIColours.MUTED_GREY


def test_input_tile_all_branches() -> None:
    """Test all controls and properties of InputTile."""
    mock_on_change = MagicMock()
    settings = GUISettings()
    settings_map: dict[str, Any] = {}
    tile = InputTile(
        settings=settings,
        settings_map=settings_map,
        on_change_handler=mock_on_change,
        header_size=18,
        font_size_default=14,
    )

    assert tile.title.value == "Input"
    assert tile.set_primer_info_panel_position is not None
    assert tile.fixed_height_primer_info_checkbox is not None
    assert tile.set_show_primer_temperature is not None
    assert tile.auto_activate_new_valid_primer is not None
    assert tile.set_tm_colour_scheme is not None
    assert tile.set_tm_dna_conc is not None
    assert tile.set_tm_dnap_conc is not None
    assert tile.set_tm_mono_salt is not None
    assert tile.set_tm_div_salt is not None
    assert tile.set_tm_dNTP_conc is not None
    assert tile.set_tm_method is not None
    assert "primer_info_panel_position" in settings_map
    assert "tm_dna_conc" in settings_map


def test_designer_1d_tile_settings() -> None:
    """Test Designer 1D tile settings in SettingsView."""
    mock_page = MagicMock(spec=ft.Page)
    settings = GUISettings()
    view = SettingsView(mock_page, settings)

    assert hasattr(view, "designer_1d_tile")
    assert view.set_designer_1d_show_tm.value is True
    assert view.set_designer_1d_show_pct_at.value is False

    # Toggle Designer 1D settings
    view.set_designer_1d_show_tm.value = False
    view.set_designer_1d_show_pct_at.value = True
    view.sync_to_state()
    assert settings["designer_1d_show_tm"] is False
    assert settings["designer_1d_show_pct_at"] is True


def test_pcr_amplicon_ranking_settings() -> None:
    """Test amplicon vertical ranking dropdown in SettingsView."""
    mock_page = MagicMock(spec=ft.Page)
    settings = GUISettings()
    view = SettingsView(mock_page, settings)

    # Test amplicon vertical ranking dropdown
    assert view.set_pcr_amplicon_ranking.value == "Position, then length"
    view.set_pcr_amplicon_ranking.value = "Quality score"
    view.sync_to_state()
    assert settings["pcr_amplicon_ranking"] == "Quality score"


def test_pcr_relative_q_score_settings() -> None:
    """Test Relative Q score checkbox in SettingsView."""
    mock_page = MagicMock(spec=ft.Page)
    settings = GUISettings()
    view = SettingsView(mock_page, settings)

    assert view.set_pcr_relative_q_score.value is True
    assert view.set_pcr_relative_q_score.label == "Relative Q score"
    view.set_pcr_relative_q_score.value = False
    view.sync_to_state()
    assert settings["pcr_relative_q_score"] is False


def test_pcr_primer_label_orientation_settings() -> None:
    """Test overview map primer label orientation dropdown in SettingsView."""
    mock_page = MagicMock(spec=ft.Page)
    settings = GUISettings()
    view = SettingsView(mock_page, settings)

    overview_section = (
        view.replication_tile.controls[0].content.controls[0].content
    )
    assert overview_section.controls[0].value == "Overview Map"
    assert overview_section.controls[1].controls == [
        view.set_pcr_divider_behaviour,
        view.set_pcr_primer_label_orientation,
    ]

    assert view.set_pcr_primer_label_orientation.value == "Horizontal"
    assert [
        opt.key for opt in view.set_pcr_primer_label_orientation.options
    ] == ["Horizontal", "45°", "Vertical"]
    view.set_pcr_primer_label_orientation.value = "45°"
    view.sync_to_state()
    assert settings["pcr_primer_label_orientation"] == "45°"
    view.set_pcr_primer_label_orientation.value = "Vertical"
    view.sync_to_state()
    assert settings["pcr_primer_label_orientation"] == "Vertical"


def test_pcr_divider_behaviour_settings() -> None:
    """Test PCR horizontal divider behaviour dropdown in SettingsView."""
    mock_page = MagicMock(spec=ft.Page)
    settings = GUISettings()
    view = SettingsView(mock_page, settings)

    assert view.set_pcr_divider_behaviour.value == "autoresize_3_4"
    assert [opt.key for opt in view.set_pcr_divider_behaviour.options] == [
        "autoresize_1_4",
        "autoresize",
        "autoresize_1_2",
        "autoresize_2_3",
        "autoresize_3_4",
        "save_proportionality",
    ]
    view.set_pcr_divider_behaviour.value = "save_proportionality"
    view.sync_to_state()
    assert settings["pcr_divider_behaviour"] == "save_proportionality"


def test_weights_cells_centred_alignment() -> None:
    """Test that all weights cells in SettingsView are centred."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    settings = GUISettings()
    view = SettingsView(mock_page, settings)

    tables = [
        view.replication_tile.score_table.table,
        view.replication_tile.match_weight_table.table,
        view.replication_tile.run_weight_table.table,
        view.dimer_tile.score_table.table,
    ]
    for table in tables:
        for row in table.rows:
            for cell in row.cells:
                container = cell.content
                assert isinstance(container, ft.Container)
                assert container.alignment == ft.Alignment(0, 0)
                if isinstance(container.content, ft.TextField):
                    assert container.content.dense is True
                    assert container.content.text_align == ft.TextAlign.CENTER
                    assert (
                        container.content.text_vertical_align
                        == ft.VerticalAlignment.CENTER
                    )


def test_dimer_weights_reset_button() -> None:
    """Test that the primer dimer weights reset button restores defaults."""
    from amplifyp.settings import DEFAULT_PRIMER_DIMER_WEIGHTS

    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    mock_on_change = MagicMock()
    settings = GUISettings()
    view = SettingsView(
        page=mock_page,
        settings=settings,
        on_change=mock_on_change,
    )

    reset_button = view.dimer_tile.score_table.reset_button
    assert isinstance(reset_button, ft.TextButton)
    assert reset_button.content == "Reset"

    # Modify a primer dimer weight and sync
    view.settings_map["pd_score_A_T"].value = "99"
    view.sync_to_state()
    assert settings["pd_score_A_T"] == "99"

    # Click the reset button (detached from a live Flet Page in unit test)
    mock_on_change.reset_mock()
    with patch.object(
        type(view.dimer_tile.score_table),
        "page",
        new_callable=lambda: property(lambda _self: None),
    ):
        view.dimer_tile.score_table._on_reset_click(None)

    assert view.settings_map["pd_score_A_T"].value == str(
        DEFAULT_PRIMER_DIMER_WEIGHTS["A", "T"]
    )
    assert settings["pd_score_A_T"] == str(
        DEFAULT_PRIMER_DIMER_WEIGHTS["A", "T"]
    )
    assert mock_on_change.called


def test_general_tile_factory_reset() -> None:
    """Test Factory Reset button and confirmation dialogue in GeneralTile."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    mock_page.overlay = []
    mock_on_change = MagicMock()
    mock_on_reset = MagicMock()

    settings = GUISettings()
    view = SettingsView(
        page=mock_page,
        settings=settings,
        on_change=mock_on_change,
        on_reset=mock_on_reset,
    )

    # Modify a wide range of preferences across sections
    settings["font_family"] = "Courier New"
    settings["dark_mode"] = True
    settings["colour_deficient"] = True
    settings["primability_cutoff"] = "0.95"
    settings["template_fixed_width"] = True
    settings["template_bases_per_line"] = "80"
    settings["auto_activate_new_valid_primer"] = True
    settings["auto_reload_on_startup"] = False
    settings["designer_1d_show_tm"] = False
    view.update_ui()

    btn = view.general_tile.factory_reset_button
    assert isinstance(btn, ft.OutlinedButton)
    assert btn.content == "Reset Settings"

    # Clicking opens the confirmation dialogue
    with patch.object(settings, "save_to_local") as mock_save:
        view.general_tile._on_factory_reset_click(MagicMock())
        dialog = view.general_tile._factory_reset_dialogue
        assert dialog is not None
        assert dialog.open is True
        assert dialog in mock_page.overlay

        # Cancel dismisses without resetting
        view.general_tile._dismiss_factory_reset()
        assert dialog.open is False
        assert settings["font_family"] == "Courier New"
        assert not mock_save.called

        # Re-open and confirm resets all preferences to factory defaults
        view.general_tile._on_factory_reset_click(MagicMock())
        assert dialog.open is True
        with patch.object(
            view.general_tile.notification_helper, "show_message"
        ) as mock_msg:
            view.general_tile._confirm_factory_reset()

        assert dialog.open is False
        assert mock_save.called
        assert mock_on_reset.called
        assert settings["font_family"] == "Roboto Mono"
        assert settings["dark_mode"] == "system"
        assert settings["colour_deficient"] is False
        assert settings["primability_cutoff"] == "0.8"
        assert settings["template_fixed_width"] is False
        assert settings["template_bases_per_line"] == "Auto"
        assert settings["auto_activate_new_valid_primer"] is False
        assert settings["auto_reload_on_startup"] is True
        assert settings["designer_1d_show_tm"] is True
        assert view.set_font_family.value == "Roboto Mono"
        assert view.set_colour_scheme.value == "System"
        assert view.set_primability_cutoff.value == "0.8"
        mock_msg.assert_called_once_with(
            "All preferences reset to their default values."
        )


def test_general_tile_factory_reset_standalone() -> None:
    """Test GeneralTile.perform_factory_reset without on_reset_callback."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    mock_page.overlay = []
    mock_on_change = MagicMock()
    mock_sync = MagicMock()
    mock_update_ui = MagicMock()

    settings = GUISettings()
    settings["font_family"] = "Consolas"
    tile = GeneralTile(
        page=mock_page,
        settings=settings,
        settings_map={},
        on_change_handler=mock_on_change,
        header_size=18,
        font_size_default=14,
        sync_to_state_callback=mock_sync,
        update_ui_callback=mock_update_ui,
    )

    with (
        patch.object(settings, "save_to_local") as mock_save,
        patch.object(tile.notification_helper, "show_message") as mock_msg,
    ):
        tile.perform_factory_reset()

    assert settings["font_family"] == "Roboto Mono"
    mock_save.assert_called_once_with(mock_page)
    mock_update_ui.assert_called_once()
    mock_on_change.assert_called_once_with(None)
    mock_msg.assert_called_once_with(
        "All preferences reset to their default values."
    )


def test_lengthwise_weight_tables_add_remove_reset_and_persistence() -> None:
    """Test LengthwiseWeightTable add, remove, reset, and GUISettings sync."""
    from amplifyp.settings import (
        DEFAULT_BASE_PAIR_WEIGHTS,
        DEFAULT_MATCH_WEIGHTS,
    )

    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    mock_on_change = MagicMock()
    settings = GUISettings()
    view = SettingsView(
        page=mock_page,
        settings=settings,
        on_change=mock_on_change,
    )

    match_tbl = view.replication_tile.match_weight_table
    # Row labels are 1-based while internal keys remain 0-based
    first_label = match_tbl.table.rows[0].cells[0].content.content.value
    assert first_label == "1"
    assert match_tbl.indices == sorted(DEFAULT_MATCH_WEIGHTS.overrides.keys())

    with patch.object(
        type(match_tbl),
        "page",
        new_callable=lambda: property(lambda _self: None),
    ):
        # 1. Add a new entry (index 15 -> displayed as "16")
        match_tbl._on_add_click(None)
        assert 15 in match_tbl.indices
        assert "match_weight_15" in view.settings_map
        assert settings["match_weight_15"] == str(
            DEFAULT_MATCH_WEIGHTS.default_weight
        )
        assert match_tbl.table.rows[-2].cells[0].content.content.value == "16"

        # 2. Remove entries down to empty (Default row remains)
        while match_tbl.indices:
            match_tbl._on_remove_click(None)
        assert match_tbl.indices == []
        assert len(match_tbl.table.rows) == 1
        assert "match_weight_0" not in settings
        assert "match_weight_15" not in settings
        assert view.get_replication_settings().match_weight.overrides == {}

        # Removing when already empty is a no-op
        match_tbl._on_remove_click(None)
        assert match_tbl.indices == []

        # 3. Round-trip added/removed indices through to_dict / from_dict
        match_tbl._on_add_click(None)  # adds index 0
        view.settings_map["match_weight_0"].value = "42"
        view.sync_to_state()
        saved = settings.to_dict()

        restored = GUISettings()
        restored.from_dict(saved)
        assert "match_weight_0" in restored
        assert "match_weight_14" not in restored
        assert restored.get_replication_settings().match_weight[0] == 42.0

        # 4. update_ui rebuilds table when settings change externally
        view.set_state(restored.to_dict())
        assert match_tbl.indices == [0]
        assert len(match_tbl.table.rows) == 2

        # 5. Per-table reset restores defaults
        match_tbl._on_reset_click(None)
        assert match_tbl.indices == sorted(
            DEFAULT_MATCH_WEIGHTS.overrides.keys()
        )
        assert settings["match_weight_0"] == str(
            DEFAULT_MATCH_WEIGHTS.overrides[0]
        )

        # 6. Base Pair Weights table also has a working Reset button
        bp_tbl = view.replication_tile.score_table
        assert bp_tbl.reset_button is not None
        view.settings_map["bp_score_G_G"].value = "123"
        view.sync_to_state()
        assert settings["bp_score_G_G"] == "123"
        with patch.object(
            type(bp_tbl),
            "page",
            new_callable=lambda: property(lambda _self: None),
        ):
            bp_tbl._on_reset_click(None)
        assert settings["bp_score_G_G"] == str(
            DEFAULT_BASE_PAIR_WEIGHTS["G", "G"]
        )


def test_gui_settings_constructor_web_get_keys_and_yaml_load() -> None:
    """Test GUISettings constructor, web get_keys fast path, and YAML load."""
    from amplifyp.gui.colours import GUIColours

    # 1. Constructor delegates to from_dict (prunes omitted lengthwise weights
    # and coerces string booleans accurately)
    init_settings = GUISettings(
        {
            "match_weight_0": "77",
            "dark_mode": "system",
            "colour_deficient": "false",
            "auto_reload_on_startup": "0",
            "amp4_compat": "no",
            "custom_extra": "ok",
        }
    )
    assert init_settings["match_weight_0"] == "77"
    assert "match_weight_14" not in init_settings
    assert init_settings["colour_deficient"] is False
    assert init_settings["auto_reload_on_startup"] is False
    assert init_settings["amp4_compat"] is False
    assert init_settings["custom_extra"] == "ok"

    # 2. Web load_from_local_async uses async get_keys fast path
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = True
    mock_prefs = MagicMock()
    mock_prefs.get_keys = AsyncMock(
        return_value=[
            "amplifyp.settings.font_family",
            "amplifyp.settings.match_weight_default",
            "amplifyp.settings.match_weight_0",
        ]
    )
    store = {
        "amplifyp.settings.font_family": "Consolas",
        "amplifyp.settings.match_weight_default": "1",
        "amplifyp.settings.match_weight_0": "88",
    }
    mock_prefs.get = AsyncMock(side_effect=lambda k: store[k])
    mock_prefs.contains_key = AsyncMock(
        side_effect=AssertionError("Should use get_keys fast path")
    )

    web_settings = GUISettings()
    with patch("flet.SharedPreferences", return_value=mock_prefs):
        asyncio.run(web_settings.load_from_local_async(mock_page))
    assert web_settings["font_family"] == "Consolas"
    assert web_settings["match_weight_0"] == "88"
    assert "match_weight_14" not in web_settings

    # 3. Loading YAML via GeneralTile inside SettingsView does not clobber
    # non-integer-MB log_rotation_max_bytes or save twice
    desktop_page = MagicMock(spec=ft.Page)
    desktop_page.web = False
    desktop_page.overlay = []
    settings = GUISettings()
    on_change_calls: list[Any] = []
    view = SettingsView(
        page=desktop_page,
        settings=settings,
        on_change=lambda e: on_change_calls.append(e),
        on_reset=lambda e: setattr(GUIColours, "dark_mode", False),
    )

    yaml_payload = yaml.safe_dump(
        {"font_family": "Courier New", "log_rotation_max_bytes": 2621440}
    )
    with (
        patch(
            "amplifyp.gui.views.settings.general_tile.pick_and_read_file",
            new=AsyncMock(return_value=yaml_payload),
        ),
        patch.object(settings, "save_to_local") as mock_save,
    ):
        asyncio.run(view.general_tile._load_settings_async(MagicMock()))

    assert settings["font_family"] == "Courier New"
    assert settings["log_rotation_max_bytes"] == 2621440
    assert len(on_change_calls) == 1
    mock_save.assert_not_called()  # delegated to controller on_change

    # 4. Resetting from Dark mode runs on_reset before update_ui and saves once
    settings["dark_mode"] = True
    GUIColours.dark_mode = True
    shared_cb = MagicMock(
        side_effect=lambda _e: setattr(GUIColours, "dark_mode", False)
    )
    view_reset = SettingsView(
        page=desktop_page,
        settings=settings,
        on_change=shared_cb,
        on_reset=shared_cb,
    )
    with patch.object(settings, "save_to_local") as mock_save_reset:
        view_reset._on_reset_handler(None)
    assert GUIColours.dark_mode is False
    assert shared_cb.call_count == 1
    assert mock_save_reset.call_count == 0

    # 5. Loading YAML with dark_mode="system" from Dark mode runs on_change
    # (apply_theme) before refreshing UI border colours
    settings["dark_mode"] = True
    GUIColours.dark_mode = True
    call_order: list[str] = []

    def _on_theme_change(_e: Any) -> None:
        call_order.append("on_change")
        GUIColours.dark_mode = False

    view_yaml_theme = SettingsView(
        page=desktop_page,
        settings=settings,
        on_change=_on_theme_change,
    )
    orig_update_ui = view_yaml_theme.update_ui

    def _tracked_update_ui() -> None:
        call_order.append(f"update_ui:dark={GUIColours.dark_mode}")
        orig_update_ui()

    view_yaml_theme.update_ui = _tracked_update_ui  # type: ignore[method-assign]
    with patch(
        "amplifyp.gui.views.settings.general_tile.pick_and_read_file",
        new=AsyncMock(return_value=yaml.safe_dump({"dark_mode": "system"})),
    ):
        asyncio.run(
            view_yaml_theme.general_tile._load_settings_async(MagicMock())
        )
    assert call_order == ["on_change", "update_ui:dark=False"]


def test_weight_tables_sync_no_defaults_and_mounted_update() -> None:
    """Test LengthwiseWeightTable sync, ScoreTable None defaults, and page."""
    from amplifyp.gui.views.settings.score_table import ScoreTable

    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    mock_on_change = MagicMock()
    settings = GUISettings()
    view = SettingsView(
        page=mock_page,
        settings=settings,
        on_change=mock_on_change,
    )

    match_tbl = view.replication_tile.match_weight_table

    # 1. sync_from_settings when indices differ rebuilds rows
    alt_settings = GUISettings(
        {
            "match_weight_0": "11",
            "match_weight_1": "22",
            "match_weight_default": "33",
        }
    )
    match_tbl.sync_from_settings(alt_settings)
    assert match_tbl.indices == [0, 1]
    assert view.settings_map["match_weight_0"].value == "11"
    assert view.settings_map["match_weight_1"].value == "22"
    assert view.settings_map["match_weight_default"].value == "33"

    # 2. sync_from_settings when indices match updates differing field values
    view.settings_map["match_weight_0"].value = "999"
    view.settings_map["match_weight_default"].value = "888"
    match_tbl.sync_from_settings(alt_settings)
    assert view.settings_map["match_weight_0"].value == "11"
    assert view.settings_map["match_weight_default"].value == "33"

    # 3. _notify_changed calls self.update() when mounted on a page
    with (
        patch.object(
            type(match_tbl),
            "page",
            new_callable=lambda: property(lambda _s: mock_page),
        ),
        patch.object(match_tbl, "update") as mock_update,
    ):
        match_tbl._notify_changed()
        mock_update.assert_called_once()

    # 4. ScoreTable._on_reset_click is a no-op when default_weights is None
    custom_map: dict[str, Any] = {"custom_A_T": ft.TextField(value="5")}
    no_default_cb = MagicMock()
    score_tbl = ScoreTable(
        label="Custom",
        row_headers=["A"],
        col_headers=["T"],
        row_label="R",
        col_label="C",
        prefix="custom",
        settings_map=custom_map,
        font_size_default=14,
        font_size_micro=10,
        font_size_table_header=14,
        default_weights=None,
        on_change_handler=no_default_cb,
    )
    score_tbl._on_reset_click(None)
    no_default_cb.assert_not_called()


def test_gui_settings_invalid_tm_method_falls_back() -> None:
    """Invalid tm_method values fall back to default TmMethod."""
    import pytest

    from amplifyp.dna import Primer

    settings = GUISettings()
    settings._settings["tm_method"] = "InvalidMethod"
    primer = Primer("ATGCGTACGT")
    assert settings.calculate_primer_tm(primer) == pytest.approx(
        GUISettings().calculate_primer_tm(primer)
    )


def test_weight_table_cell_text_style_refreshes_on_theme_change() -> None:
    """New/reset/synced weight cells adopt updated theme colour and size."""
    from amplifyp.gui.colours import GUIColours

    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    settings = GUISettings()
    view = SettingsView(
        page=mock_page,
        settings=settings,
        on_change=MagicMock(),
    )
    match_tbl = view.replication_tile.match_weight_table
    assert (
        view.settings_map["match_weight_0"].text_style.color
        == GUIColours.DIAGRAM_BLACK
    )

    custom_colour = GUIColours.PRIMARY
    assert custom_colour != GUIColours.DIAGRAM_BLACK

    with patch.object(
        type(GUIColours),
        "DIAGRAM_BLACK",
        new_callable=lambda: property(lambda _cls: custom_colour),
    ):
        # update_colours refreshes the cached TextStyle
        match_tbl.update_colours(match_tbl.table)
        assert match_tbl._cell_text_style.color == custom_colour

        # Newly added cell uses refreshed text colour
        with patch.object(
            type(match_tbl),
            "page",
            new_callable=lambda: property(lambda _s: None),
        ):
            match_tbl._on_add_click(None)
        added_key = f"match_weight_{match_tbl.indices[-1]}"
        assert view.settings_map[added_key].text_style.color == custom_colour

        # Reset cells use refreshed text colour
        with patch.object(
            type(match_tbl),
            "page",
            new_callable=lambda: property(lambda _s: None),
        ):
            match_tbl._on_reset_click(None)
        assert (
            view.settings_map["match_weight_0"].text_style.color
            == custom_colour
        )

    # Rebuilt cells via sync_from_settings pick up restored colour and new size
    match_tbl.font_size_default = 18
    alt_settings = GUISettings(
        {"match_weight_0": "10", "match_weight_default": "1"}
    )
    match_tbl.sync_from_settings(alt_settings)
    rebuilt_style = view.settings_map["match_weight_0"].text_style
    assert rebuilt_style.color == GUIColours.DIAGRAM_BLACK
    assert rebuilt_style.size == 18
