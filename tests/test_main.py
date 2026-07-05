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
#

"""Tests for the GUI app."""

from unittest.mock import MagicMock

import flet
import flet as ft
import pytest

from amplifyp.gui.app import main as gui_main
from main import main


def test_gui_main() -> None:
    """Test the main function of the GUI app."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.window = MagicMock()

    main(mock_page)

    assert mock_page.title == "AmplifyP"
    assert mock_page.vertical_alignment == ft.MainAxisAlignment.START
    mock_page.add.assert_called()


def test_gui_app_main() -> None:
    """Test the main function imported directly from amplifyp.gui.app."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.window = MagicMock()

    gui_main(mock_page)

    assert mock_page.title == "AmplifyP"
    assert mock_page.vertical_alignment == ft.MainAxisAlignment.START
    mock_page.add.assert_called()


def test_window_close_warning_dialog() -> None:
    """Test setup, trigger, and dismiss of close warning dialog.

    Simulates a hot reload scenario as well.
    """
    mock_page = MagicMock(spec=ft.Page)
    mock_page.window = MagicMock()
    mock_page.overlay = []
    mock_page.web = False

    # Call main first time (initial startup)
    gui_main(mock_page)
    assert mock_page.window.prevent_close is True

    # Call main second time (simulating a hot reload on the same Page instance)
    gui_main(mock_page)

    # 1. Verify prevent_close is still enabled
    assert mock_page.window.prevent_close is True

    # 2. Verify event handler is registered
    assert mock_page.window.on_event is not None

    # 3. Simulate close event
    mock_event = MagicMock()
    mock_event.data = "close"
    mock_page.window.on_event(mock_event)

    # Verify dialog is created, appended to overlay and opened
    dialog = next(
        (c for c in mock_page.overlay if isinstance(c, ft.AlertDialog)), None
    )
    assert dialog is not None
    assert dialog.open is True
    mock_page.update.assert_called()

    # Reset update mock for dismissal check
    mock_page.update.reset_mock()

    # 4. Simulate clicking "No" (confirm_dismiss)
    # The actions of AlertDialog are: Yes (index 0) and No (index 1)
    no_button = dialog.actions[1]
    no_button.on_click(MagicMock())

    # Verify dialog is closed
    assert dialog.open is False
    mock_page.update.assert_called_once()


def test_cli_default(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test cli() with default desktop settings."""
    import main as main_module

    mock_run = MagicMock()
    monkeypatch.setattr(flet, "run", mock_run)

    main_module.cli([])

    assert main_module.state_file is None
    assert main_module.auto_close is False
    mock_run.assert_called_once()
    kwargs = mock_run.call_args.kwargs
    assert kwargs["view"] == flet.AppView.FLET_APP
    assert kwargs["port"] == 0


def test_cli_web(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test cli() with --web flag."""
    import main as main_module

    mock_run = MagicMock()
    monkeypatch.setattr(flet, "run", mock_run)

    main_module.cli(["--web"])

    assert main_module.state_file is None
    assert main_module.auto_close is False
    mock_run.assert_called_once()
    kwargs = mock_run.call_args.kwargs
    assert kwargs["view"] == flet.AppView.WEB_BROWSER
    assert kwargs["port"] == 34521
