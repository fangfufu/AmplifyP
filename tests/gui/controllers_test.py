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

"""Tests for GUI controllers and app lifecycle."""

from unittest.mock import AsyncMock, MagicMock, patch

import flet as ft
import pytest
from flet.messaging.session import Session

from amplifyp.gui.app import _patch_flet_session
from amplifyp.gui.app import main as app_main
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.controller import GUIController
from amplifyp.gui.controllers.navigation import NavigationManager
from amplifyp.gui.controllers.theme import ThemeManager
from amplifyp.gui.controllers.updater import UpdateManager
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.views.settings.primer_list_tile import PrimerListTile


def test_navigation_manager_flow() -> None:
    """Test NavigationManager setup, view switching, and action handlers."""
    mock_ctrl = MagicMock()
    mock_ctrl.settings = GUISettings()
    mock_ctrl.page = MagicMock(spec=ft.Page)
    mock_ctrl.page.controls = []
    mock_ctrl.pcr_button_ref = ft.Ref[ft.FilledButton]()
    mock_ctrl.dimers_button_ref = ft.Ref[ft.FilledButton]()
    mock_ctrl.input_view = MagicMock()
    mock_ctrl.settings_view = MagicMock()
    mock_ctrl.about_view = MagicMock()
    mock_ctrl.designer_view = MagicMock()
    mock_ctrl.designer_2d_view = MagicMock()
    mock_ctrl.pcr_view = MagicMock()
    mock_ctrl.dimers_view = MagicMock()
    mock_ctrl.view_container = MagicMock()
    mock_ctrl.input_view_dirty = True

    nav_manager = NavigationManager(mock_ctrl)
    nav_manager.setup_navigation_controls()

    assert mock_ctrl.header is not None
    assert mock_ctrl.header_container is not None

    # Switch to dirty input_view
    nav_manager.switch_view(MagicMock(), mock_ctrl.input_view)
    assert mock_ctrl.input_view_dirty is False
    mock_ctrl.input_view.update_ui.assert_called()
    assert mock_ctrl.page.on_resize == mock_ctrl.input_view._handle_resize

    # Switch to pcr_view
    nav_manager.switch_view(MagicMock(), mock_ctrl.pcr_view)
    assert mock_ctrl.page.on_resize == mock_ctrl.pcr_view._handle_resize

    # Switch to settings_view
    nav_manager.switch_view(MagicMock(), mock_ctrl.settings_view)
    assert mock_ctrl.page.on_resize is None

    # on_pcr_click when run_pcr returns False -> reverts to input_view
    mock_ctrl.pcr_view.run_pcr.return_value = False
    nav_manager.on_pcr_click(MagicMock())
    assert mock_ctrl.view_container.content == mock_ctrl.input_view

    # on_pcr_click when run_pcr returns True -> stays on pcr_view
    mock_ctrl.pcr_view.run_pcr.return_value = True
    nav_manager.on_pcr_click(MagicMock())
    assert mock_ctrl.view_container.content == mock_ctrl.pcr_view

    # on_dimers_click when run_analysis returns False -> reverts to input_view
    mock_ctrl.dimers_view.run_analysis.return_value = False
    nav_manager.on_dimers_click(MagicMock())
    assert mock_ctrl.view_container.content == mock_ctrl.input_view

    # on_dimers_click when run_analysis returns True -> stays on dimers_view
    mock_ctrl.dimers_view.run_analysis.return_value = True
    nav_manager.on_dimers_click(MagicMock())
    assert mock_ctrl.view_container.content == mock_ctrl.dimers_view


def test_theme_manager_and_brightness_changes() -> None:
    """Test ThemeManager theme mode resolution and brightness changes."""

    mock_ctrl = MagicMock()
    mock_ctrl.settings = GUISettings()
    mock_ctrl.page = MagicMock(spec=ft.Page)
    mock_ctrl.header_container = MagicMock()
    mock_ctrl.view_container = MagicMock()
    mock_ctrl.input_view = MagicMock()
    mock_ctrl.settings_view = MagicMock()
    mock_ctrl.pcr_view = MagicMock()
    mock_ctrl.dimers_view = MagicMock()

    theme_manager = ThemeManager(mock_ctrl)

    # 1. System dark mode with dark platform brightness
    mock_ctrl.settings["dark_mode"] = "system"
    mock_ctrl.page.platform_brightness = "dark"
    theme_manager.apply_theme()
    assert mock_ctrl.page.theme_mode == ft.ThemeMode.SYSTEM
    assert GUIColours.dark_mode is True

    # 2. System dark mode with light platform brightness
    mock_ctrl.page.platform_brightness = "light"
    theme_manager.apply_theme()
    assert GUIColours.dark_mode is False

    # 3. Explicit dark mode
    mock_ctrl.settings["dark_mode"] = True
    theme_manager.apply_theme()
    assert mock_ctrl.page.theme_mode == ft.ThemeMode.DARK
    assert GUIColours.dark_mode is True

    # 4. Explicit light mode
    mock_ctrl.settings["dark_mode"] = False
    theme_manager.apply_theme()
    assert mock_ctrl.page.theme_mode == ft.ThemeMode.LIGHT
    assert GUIColours.dark_mode is False

    # 5. on_platform_brightness_change across various views
    mock_ctrl.view_container.content = mock_ctrl.input_view
    theme_manager.on_platform_brightness_change()
    mock_ctrl.input_view.update_ui.assert_called()

    mock_ctrl.view_container.content = mock_ctrl.settings_view
    theme_manager.on_platform_brightness_change()
    assert mock_ctrl.input_view_dirty is True
    mock_ctrl.settings_view.update_ui.assert_called()

    mock_ctrl.view_container.content = mock_ctrl.pcr_view
    theme_manager.on_platform_brightness_change()
    mock_ctrl.pcr_view.run_pcr.assert_called_with(keep_cards=True)

    mock_ctrl.view_container.content = mock_ctrl.dimers_view
    theme_manager.on_platform_brightness_change()
    mock_ctrl.dimers_view.run_analysis.assert_called()


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_update_manager_async() -> None:
    """Test UpdateManager version check scheduling and update notification."""
    mock_ctrl = MagicMock()
    mock_ctrl.settings = GUISettings()
    mock_ctrl.header = MagicMock()
    mock_ctrl.page = MagicMock()
    update_manager = UpdateManager(mock_ctrl)

    # on_update_found
    update_manager.on_update_found("v2.1.0")
    mock_ctrl.header.set_update_available.assert_called_with("v2.1.0")

    # check_updates_async skipped if frequency is Never
    mock_ctrl.settings["version_checking_frequency"] = "Never"
    await update_manager.check_updates_async()

    # check_updates_async with update found
    mock_ctrl.settings["version_checking_frequency"] = "Once per Month"
    mock_ctrl.settings["last_version_check_timestamp"] = "invalid_timestamp"
    with patch(
        "amplifyp.gui.utils.system.fetch_latest_release_version",
        return_value="v99.0.0",
    ):
        await update_manager.check_updates_async()
        mock_ctrl.header.set_update_available.assert_called_with("v99.0.0")


def test_flet_session_patch_and_edge_cases() -> None:
    """Test _patch_flet_session idempotency and message handling behavior."""
    # Ensure patched once
    _patch_flet_session()

    # Create mock session simulating Flet internal state
    session = MagicMock(spec=Session)
    setattr(session, f"_{Session.__name__}__index", {10: "control_10"})
    evt_registered = MagicMock()
    evt_unregistered = MagicMock()
    method_calls = {"call_1": evt_registered, "call_2": evt_unregistered}
    setattr(session, f"_{Session.__name__}__method_calls", method_calls)
    setattr(session, f"_{Session.__name__}__method_call_results", {})

    # 1. Invoke method result for registered control
    Session.handle_invoke_method_results(  # pyright: ignore[reportAttributeAccessIssue]
        session, 10, "call_1", "result_ok", None
    )
    evt_registered.set.assert_called_once()

    # 2. Invoke method result for unregistered control
    Session.handle_invoke_method_results(  # pyright: ignore[reportAttributeAccessIssue]
        session, 999, "call_2", "result_orphan", None
    )
    evt_unregistered.set.assert_called_once()

    # 3. Missing call_id returns safely without error
    Session.handle_invoke_method_results(  # pyright: ignore[reportAttributeAccessIssue]
        session, 10, "unknown_call", None, None
    )
    Session.handle_invoke_method_results(  # pyright: ignore[reportAttributeAccessIssue]
        session, 999, "unknown_call", None, None
    )

    # 4. Exception in _patch_flet_session gracefully handled
    Session._amplifyp_patched = False
    with patch(
        "amplifyp.gui.app.getattr",
        side_effect=Exception("Failed to patch"),
    ):
        _patch_flet_session()


def test_app_main_options() -> None:
    """Test app.py main entry point with custom window dimensions."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    mock_page.window = MagicMock()
    mock_page.controls = []

    app_main(
        mock_page,
        window_width=1400,
        window_height=900,
        auto_close=True,
    )

    assert mock_page.window.width == 1400
    assert mock_page.window.height == 900
    assert mock_page.window.prevent_close is False


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_gui_controller_methods_and_delegations() -> None:
    """Test GUIController initialisation and delegation methods."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = True
    mock_page.run_javascript = MagicMock()
    mock_page.window = MagicMock()
    mock_page.controls = []

    controller = GUIController(
        page=mock_page,
        auto_close=False,
        state_file="custom_state.yaml",
    )
    controller.initialise()

    # Trigger input change and stop editing callbacks
    if (change_cb := controller.input_view.on_change) is not None:
        change_cb(None)
    cb = controller.input_view.on_stop_editing_callback

    if cb is not None:
        with patch.object(controller, "save_last_state") as mock_save:
            cb(None)
            mock_save.assert_called_once()

    # Desktop configuration
    mock_page.web = False
    controller._configure_page_and_window()

    # View switching and callbacks
    controller.on_platform_brightness_change(None)
    controller.switch_view(MagicMock(), controller.settings_view)
    controller.on_pcr_click(MagicMock())
    controller.on_dimers_click(MagicMock())

    with patch.object(controller._update_manager, "on_update_found"):
        controller.on_update_found("v1.0.0")

    with patch.object(
        controller._update_manager,
        "check_updates_async",
        new_callable=AsyncMock,
    ):
        await controller.check_updates_async()

    # Settings change handling across views
    controller.view_container.content = controller.input_view
    controller.on_settings_change()

    controller.view_container.content = controller.settings_view
    controller.on_settings_change()

    controller.view_container.content = controller.pcr_view
    controller.on_settings_change()

    controller.view_container.content = controller.dimers_view
    controller.on_settings_change()

    # Primer list tile reposition branch in settings change
    mock_ev = MagicMock()
    controller.settings_view.primer_list_tile = MagicMock(spec=PrimerListTile)
    pos_ctrl = MagicMock()
    set_pos = controller.settings_view.primer_list_tile
    set_pos.set_primer_info_panel_position = pos_ctrl
    mock_ev.control = pos_ctrl

    controller.view_container.content = controller.input_view
    controller.on_settings_change(mock_ev)

    # Exception safety branch when checking primer list tile
    class BadSettingsView:
        @property
        def primer_list_tile(self) -> None:
            raise RuntimeError("tile lookup failed")

    controller.settings_view = BadSettingsView()  # type: ignore[assignment]
    controller.on_settings_change(mock_ev)

    # PCR button state updates with validation errors and active primers
    controller.input_data.template = "ATGCATGC"
    controller.input_data.primers = [
        {"name": "P1", "seq": "ATGC", "active": True}
    ]
    controller.input_view.primer_input.validation_errors = [
        {"name": None, "seq": None}
    ]
    controller.update_pcr_button_state(sync=False)

    # With validation error on active primer
    controller.input_view.primer_input.validation_errors = [
        {"name": "Invalid primer", "seq": None}
    ]
    controller.update_pcr_button_state(sync=False)

    # Submodule helper delegations
    with (
        patch("amplifyp.gui.utils.data_helpers.get_last_state_path"),
        patch("amplifyp.gui.utils.data_helpers.save_last_state"),
        patch("amplifyp.gui.utils.data_helpers.load_last_state"),
        patch("amplifyp.gui.utils.data_helpers.restore_state_from_file"),
        patch("amplifyp.gui.utils.data_helpers.apply_parsed_state"),
        patch(
            "amplifyp.gui.utils.data_helpers.save_state",
            new_callable=AsyncMock,
        ),
        patch(
            "amplifyp.gui.utils.data_helpers.load_state",
            new_callable=AsyncMock,
        ),
        patch("amplifyp.gui.utils.data_helpers.confirm_clear"),
        patch("amplifyp.gui.utils.data_helpers.dismiss_clear"),
        patch("amplifyp.gui.utils.data_helpers.clear_all"),
        patch("amplifyp.gui.utils.system.confirm_dismiss"),
        patch(
            "amplifyp.gui.utils.system.confirm_exit_async",
            new_callable=AsyncMock,
        ),
        patch("amplifyp.gui.utils.system.confirm_exit"),
        patch(
            "amplifyp.gui.utils.system.restore_state_and_auto_close_async",
            new_callable=AsyncMock,
        ),
        patch(
            "amplifyp.gui.utils.system.auto_close_and_quit_delayed",
            new_callable=AsyncMock,
        ),
        patch("amplifyp.gui.utils.system.on_window_event"),
        patch("amplifyp.gui.utils.gui_helpers.handle_keyboard_event"),
    ):
        controller._get_last_state_path()
        controller.save_last_state()
        controller.load_last_state()
        await controller._load_last_state_async()
        controller._restore_state_from_file("dummy.yaml")
        controller._apply_parsed_state({})
        await controller.save_state(MagicMock())
        await controller.load_state(MagicMock())
        controller.confirm_dismiss(MagicMock())
        await controller.confirm_exit_async()
        controller.confirm_exit(MagicMock())
        await controller._restore_state_and_auto_close_async()
        await controller._auto_close_and_quit_delayed(MagicMock())
        controller.on_window_event(MagicMock())
        controller._confirm_clear(MagicMock())
        controller._dismiss_clear(MagicMock())
        controller.clear_all(MagicMock())
        controller._on_keyboard_event(MagicMock())
