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

    # 1. Properties
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
    assert view.set_colour_deficient is not None
    assert view.set_improved_visualisation is not None
    assert view.set_show_primer_temperature is not None
    assert view.set_tm_colour_scheme is not None
    assert view.set_designer_2d_colour_scheme is not None
    assert view.set_designer_2d_show_rev_fwd is not None
    assert view.set_version_checking_frequency is not None
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
    view._on_reset_handler(mock_event)
    assert mock_on_reset.called

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
