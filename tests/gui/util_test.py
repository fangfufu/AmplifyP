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

"""Tests for GUI utility functions."""

import subprocess
import sys
from typing import Any, cast
from unittest.mock import MagicMock

import flet as ft
import pytest

from amplifyp.gui.utils.data_helpers import clean_sequence, format_sequence
from amplifyp.gui.utils.gui_helpers import show_error_dialog


def test_clean_sequence() -> None:
    """Test clean_sequence removes whitespaces and newlines."""
    assert clean_sequence(" A T \n G C ") == "ATGC"
    assert clean_sequence("") == ""


def test_format_sequence() -> None:
    """Test format_sequence wraps long sequences at wrap_length."""
    seq = "ATGC" * 10  # 40 bp
    formatted = format_sequence(seq, wrap_length=10)
    assert formatted == "ATGCATGCAT\nGCATGCATGC\nATGCATGCAT\nGCATGCATGC"


def test_show_error_dialog() -> None:
    """Test show_error_dialog correctly manages overlay and OK click."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.overlay = []

    show_error_dialog(mock_page, "Error Title", "Test Message")

    # 1. Verify dialog was added to overlay and set to open
    assert len(mock_page.overlay) == 1
    dialog = mock_page.overlay[0]
    assert isinstance(dialog, ft.AlertDialog)
    assert dialog.open is True

    # 2. Click the OK button action
    assert len(dialog.actions) == 1
    ok_button = dialog.actions[0]
    ok_button.on_click(MagicMock())

    # Verify clicking OK closes the dialog (open = False)
    assert dialog.open is False

    # 3. Simulate Flet's on_dismiss event triggered on dismissal
    dialog.on_dismiss(MagicMock())

    # Verify dialog was cleaned up and removed from overlay
    assert dialog not in mock_page.overlay


def test_get_version_and_sha() -> None:
    """Test get_git_sha, get_full_sha, and get_version."""
    from unittest.mock import MagicMock, patch

    from amplifyp.gui.utils.system import get_full_sha, get_git_sha, get_version

    # Test under mocked 'amplifyp.gui.git_sha' module
    mock_git_sha = MagicMock()
    mock_git_sha.GIT_SHA = "pkgsha"
    mock_git_sha.GIT_FULL_SHA = "pkgfullsha"
    with patch.dict(sys.modules, {"amplifyp.gui.git_sha": mock_git_sha}):
        assert get_git_sha() == "pkgsha"
        assert get_full_sha() == "pkgfullsha"

    # Test under mocked 'js' module
    mock_js = MagicMock()
    mock_js.window = MagicMock()
    mock_js.window.__APP_SHA__ = "mockedsha"

    # Hide amplifyp.gui.git_sha to test js fallback
    with patch.dict(sys.modules, {"amplifyp.gui.git_sha": None, "js": mock_js}):
        # Trigger module reload if needed, or intercept imports
        assert get_git_sha() == "mockedsha"
        assert get_full_sha() == "mockedsha"

    # Test get_version uses __version__ defined in __init__.py
    from amplifyp import __version__

    version_str = get_version()
    assert __version__ in version_str


def test_handle_keyboard_event_template_copy() -> None:
    """Test Ctrl+C / Cmd+C copies sequence without linebreaks."""
    from unittest.mock import patch

    from amplifyp.gui.utils.gui_helpers import handle_keyboard_event

    mock_controller = MagicMock()
    mock_input_view = MagicMock()
    mock_template_input = MagicMock()
    mock_template_sequence = MagicMock()

    mock_controller.input_view = mock_input_view
    mock_controller.view_container.content = mock_input_view
    mock_input_view.template_input = mock_template_input
    mock_template_input.template_sequence = mock_template_sequence
    mock_input_view._currently_focused_control = mock_template_sequence

    mock_template_sequence.value = "ATGC\nATGC\nATGC"
    mock_template_sequence.selection = None

    mock_event = MagicMock(spec=ft.KeyboardEvent)
    mock_event.key = "c"
    mock_event.ctrl = True
    mock_event.meta = False

    mock_controller.page.web = False

    with patch("pyperclip.copy") as mock_pyperclip_copy:
        handle_keyboard_event(mock_controller, mock_event)
        mock_pyperclip_copy.assert_called_once_with("ATGCATGCATGC")


def test_handle_keyboard_event_template_copy_newline_only() -> None:
    """Test Ctrl+C with newline-only template does not write to clipboard."""
    from unittest.mock import patch

    from amplifyp.gui.utils.gui_helpers import handle_keyboard_event

    mock_controller = MagicMock()
    mock_input_view = MagicMock()
    mock_template_input = MagicMock()
    mock_template_sequence = MagicMock()

    mock_controller.input_view = mock_input_view
    mock_controller.view_container.content = mock_input_view
    mock_input_view.template_input = mock_template_input
    mock_template_input.template_sequence = mock_template_sequence
    mock_input_view._currently_focused_control = mock_template_sequence

    mock_template_sequence.value = "\n\n\n"
    mock_template_sequence.selection = None

    mock_event = MagicMock(spec=ft.KeyboardEvent)
    mock_event.key = "c"
    mock_event.ctrl = True
    mock_event.meta = False

    mock_controller.page.web = False

    with patch("pyperclip.copy") as mock_pyperclip_copy:
        handle_keyboard_event(mock_controller, mock_event)
        mock_pyperclip_copy.assert_not_called()


def test_git_fallback_to_dot_git() -> None:
    """Test git utility fallback to reading .git directory directly."""
    from typing import Any
    from unittest.mock import mock_open, patch

    # Mock subprocess.run to raise SubprocessError.
    # Mock sys.modules to remove git_sha and js.
    # Mock os.path.exists and open to simulate .git.
    with (
        patch("subprocess.run") as mock_run,
        patch("os.path.exists") as mock_exists,
        patch(
            "builtins.open", mock_open(read_data="ref: refs/heads/main\n")
        ) as mock_file_open,
        patch.dict(sys.modules, {"amplifyp.gui.git_sha": None, "js": None}),
    ):
        mock_run.side_effect = subprocess.SubprocessError("git failed")

        # Make os.path.exists return True for HEAD and ref file.
        def exists_side_effect(path: str) -> bool:
            return ".git" in path or "git-sha" in path

        mock_exists.side_effect = exists_side_effect

        # Let's customize the file contents read
        file_contents = {
            "HEAD": "ref: refs/heads/test-branch\n",
            "test-branch": "a1b2c3d4e5f6g7h8i9j0a1b2c3d4e5f6g7h8i9j0\n",
        }

        # We can implement a side effect for mock_file_open
        def open_side_effect(file: Any, *args: Any, **kwargs: Any) -> Any:
            filename = str(file)
            for key, val in file_contents.items():
                if filename.endswith(key):
                    # Return mock file handle supporting context manager.
                    mo = mock_open(read_data=val)
                    return mo.return_value
            raise FileNotFoundError(f"File not found: {filename}")

        mock_file_open.side_effect = open_side_effect

        from amplifyp.gui.utils.system import get_full_sha, get_git_sha

        # Verify it resolves correctly
        assert get_git_sha() == "a1b2c3d"
        assert get_full_sha() == "a1b2c3d4e5f6g7h8i9j0a1b2c3d4e5f6g7h8i9j0"

        # Detached HEAD state
        file_contents_detached = {
            "HEAD": "9999999999999999999999999999999999999999\n"
        }

        def open_side_effect_detached(
            file: Any, *args: Any, **kwargs: Any
        ) -> Any:
            filename = str(file)
            for key, val in file_contents_detached.items():
                if filename.endswith(key):
                    mo = mock_open(read_data=val)
                    return mo.return_value
            raise FileNotFoundError(f"File not found: {filename}")

        mock_file_open.side_effect = open_side_effect_detached

        assert get_git_sha() == "9999999"
        assert get_full_sha() == "9999999999999999999999999999999999999999"


def test_serialise_state_multiline() -> None:
    """Test serialise_state uses | block style for multiline strings."""
    from amplifyp.gui.util import serialise_state

    state: dict[str, object] = {
        "single": "single_line",
        "multi": "line1\nline2\nline3\n",
    }
    yaml_out = serialise_state(state)
    assert "|" in yaml_out
    assert "line1" in yaml_out


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_data_helpers_and_system_utilities(tmp_path: Any) -> None:
    """Test data helpers, UI widgets, file ops, and system lifecycle."""
    from unittest.mock import AsyncMock, mock_open, patch

    from amplifyp.gui.controller import GUIController
    from amplifyp.gui.utils.data_helpers import (
        _resolve_font_family,
        clear_all,
        confirm_clear,
        create_overlapped_sequence_view,
        dismiss_clear,
        load_last_state,
        load_state,
        pick_and_read_file,
        restore_state_from_file,
        save_and_write_file,
        save_last_state,
        save_state,
    )
    from amplifyp.gui.utils.gui_helpers import (
        BorderedCheckbox,
        Debouncer,
        handle_keyboard_event,
        initialise_score_fields,
    )
    from amplifyp.gui.utils.system import (
        _get_sha,
        auto_close_and_quit_delayed,
        capture_all_views_async,
        capture_view_screenshot_async,
        confirm_dismiss,
        confirm_exit,
        confirm_exit_async,
        get_version,
        on_window_event,
        restore_state_and_auto_close_async,
    )
    from amplifyp.gui.views.input.primer.row import PrimerRow

    # 1. Fonts and overlapped sequence views
    assert _resolve_font_family("") == "Roboto Mono"
    assert _resolve_font_family("Consolas") == "Consolas"
    assert _resolve_font_family("roboto mono") == "Roboto Mono"
    assert _resolve_font_family("custom_font") == "custom_font"

    t_dimer = create_overlapped_sequence_view(
        "TOP",
        "MID",
        "BOT",
        is_dimer=True,
        top_name_line="TOP_NAME",
        bottom_name_line="BOT_NAME",
    )
    assert len(t_dimer.spans) == 5

    t_normal_1 = create_overlapped_sequence_view(
        "TOP", "MID", "SINGLE_BOT", is_dimer=False
    )
    assert len(t_normal_1.spans) == 5

    _ = create_overlapped_sequence_view("TOP", "MID", "B1\nB2", is_dimer=False)
    _ = create_overlapped_sequence_view(
        "TOP", "MID", "B1\nB2\nB3", is_dimer=False
    )

    # 2. Debouncer and UI widgets
    debouncer = Debouncer(delay_seconds=0.01)
    triggered = False

    def on_debounced() -> None:
        nonlocal triggered
        triggered = True

    debouncer.trigger(on_debounced)
    debouncer.cancel()

    with patch(
        "threading.Timer.start", side_effect=RuntimeError("Thread error")
    ):
        debouncer.trigger(on_debounced)
        assert triggered is True

    checkbox = BorderedCheckbox(label="Option", value=False)
    assert checkbox.value is False
    checkbox.value = "true"
    assert checkbox.value is True
    checkbox.value = True
    assert checkbox.value is True

    score_map: dict[str, Any] = {}
    initialise_score_fields(
        settings_map=score_map,
        prefix="score",
        row_headers=["A"],
        col_headers=["T"],
        on_change_handler=MagicMock(),
        font_size=12,
    )
    assert "score_A_T" in score_map

    # 3. File I/O helpers
    mock_page = MagicMock(spec=ft.Page)
    mock_page.services = []
    mock_page.overlay = []
    mock_page.web = False
    mock_page.window = MagicMock()
    mock_page.controls = []
    mock_page.client_storage = MagicMock()
    mock_page.client_storage.contains_key.return_value = True
    mock_page.client_storage.get.return_value = {
        "template": "ATGC",
        "primers": [],
    }

    notifier = MagicMock()

    # pick_and_read_file variations
    with patch(
        "flet.FilePicker.pick_files",
        new_callable=AsyncMock,
        return_value=None,
    ):
        assert (
            await pick_and_read_file(mock_page, "Title", ["yaml"], notifier)
            is None
        )

    mock_file_bytes = MagicMock()
    mock_file_bytes.bytes = b"hello bytes"
    with patch(
        "flet.FilePicker.pick_files",
        new_callable=AsyncMock,
        return_value=[mock_file_bytes],
    ):
        assert (
            await pick_and_read_file(mock_page, "Title", ["yaml"], notifier)
            == "hello bytes"
        )

    mock_file_path = MagicMock()
    mock_file_path.bytes = None
    mock_file_path.path = "/some/file.txt"
    with (
        patch(
            "flet.FilePicker.pick_files",
            new_callable=AsyncMock,
            return_value=[mock_file_path],
        ),
        patch("builtins.open", mock_open(read_data="file text")),
    ):
        assert (
            await pick_and_read_file(mock_page, "Title", ["yaml"], notifier)
            == "file text"
        )

    mock_file_none = MagicMock()
    mock_file_none.bytes = None
    mock_file_none.path = None
    with patch(
        "flet.FilePicker.pick_files",
        new_callable=AsyncMock,
        return_value=[mock_file_none],
    ):
        assert (
            await pick_and_read_file(mock_page, "Title", ["yaml"], notifier)
            is None
        )

    with patch(
        "flet.FilePicker.pick_files",
        new_callable=AsyncMock,
        side_effect=OSError("io error"),
    ):
        assert (
            await pick_and_read_file(mock_page, "Title", ["yaml"], notifier)
            is None
        )

    # save_and_write_file variations
    mock_page.web = True
    with patch(
        "flet.FilePicker.save_file", new_callable=AsyncMock, return_value=None
    ):
        assert (
            await save_and_write_file(
                mock_page, "Save", "f.txt", ["txt"], "content", notifier
            )
            is True
        )

    mock_page.web = False
    with patch(
        "flet.FilePicker.save_file", new_callable=AsyncMock, return_value=None
    ):
        assert (
            await save_and_write_file(
                mock_page, "Save", "f.txt", ["txt"], "content", notifier
            )
            is False
        )

    with (
        patch(
            "flet.FilePicker.save_file",
            new_callable=AsyncMock,
            return_value="/saved.txt",
        ),
        patch("builtins.open", mock_open()),
    ):
        assert (
            await save_and_write_file(
                mock_page, "Save", "f.txt", ["txt"], "content", notifier
            )
            is True
        )

    with patch(
        "flet.FilePicker.save_file",
        new_callable=AsyncMock,
        side_effect=OSError("write error"),
    ):
        assert (
            await save_and_write_file(
                mock_page, "Save", "f.txt", ["txt"], "content", notifier
            )
            is False
        )

    # 4. State persistence & restoration
    ctrl = GUIController(page=mock_page)
    ctrl.initialise()

    mock_page.web = True
    save_last_state(ctrl)
    load_last_state(ctrl)

    ctrl.settings["auto_reload_on_startup"] = False
    save_last_state(ctrl)
    load_last_state(ctrl)
    ctrl.settings["auto_reload_on_startup"] = True

    mock_page.web = False
    with patch("builtins.open", side_effect=OSError("open error")):
        save_last_state(ctrl)
        load_last_state(ctrl)

    with patch("builtins.open", mock_open(read_data="- list item")):
        restore_state_from_file(ctrl, "invalid.yaml")
    with patch("builtins.open", side_effect=OSError("file missing")):
        restore_state_from_file(ctrl, "missing.yaml")

    ctrl.filepicker_open = True
    await save_state(ctrl, MagicMock())
    await load_state(ctrl, MagicMock())
    ctrl.filepicker_open = False

    with patch(
        "amplifyp.gui.utils.data_helpers.save_and_write_file",
        new_callable=AsyncMock,
    ):
        await save_state(ctrl, MagicMock())

    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new_callable=AsyncMock,
        return_value="- list item",
    ):
        await load_state(ctrl, MagicMock())

    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new_callable=AsyncMock,
        return_value="input: {template: ATGC}",
    ):
        await load_state(ctrl, MagicMock())

    confirm_clear(ctrl, MagicMock())
    dismiss_clear(ctrl, MagicMock())
    clear_all(ctrl, MagicMock())

    # 5. Keyboard events and navigation
    ctrl.view_container.content = ctrl.settings_view
    handle_keyboard_event(ctrl, MagicMock())
    ctrl.view_container.content = ctrl.input_view

    # Copy template sequence
    ctrl.input_view._currently_focused_control = (
        ctrl.input_view.template_input.template_sequence
    )
    ctrl.input_view.template_input.template_sequence.value = "A T G C\n"
    ctrl.input_view.template_input.template_sequence.selection = (
        ft.TextSelection(base_offset=0, extent_offset=3)
    )

    ev_copy = MagicMock()
    ev_copy.key = "c"
    ev_copy.ctrl = True
    ev_copy.meta = False

    with patch("pyperclip.copy") as mock_clip:
        handle_keyboard_event(ctrl, ev_copy)
        mock_clip.assert_called_with("AT")

    mock_page.web = True
    mock_page.run_javascript = MagicMock()
    handle_keyboard_event(ctrl, ev_copy)
    mock_page.run_javascript.assert_called()
    mock_page.web = False

    # Primer row keyboard navigation
    row0 = MagicMock(spec=PrimerRow)
    row0.idx = 0
    row0.checkbox = MagicMock(value=True)
    row0.name_field = ft.TextField(
        value="Fwd", data={"idx": 0, "field": "name", "cursor_pos": 3}
    )
    row0.name_field.update = MagicMock()
    row0.name_field.focus = MagicMock()
    row0.seq_field = ft.TextField(
        value="ATGC", data={"idx": 0, "field": "seq", "cursor_pos": 0}
    )
    row0.seq_field.update = MagicMock()
    row0.seq_field.focus = MagicMock()

    row1 = MagicMock(spec=PrimerRow)
    row1.idx = 1
    row1.checkbox = MagicMock(value=True)
    row1.name_field = ft.TextField(
        value="Rev", data={"idx": 1, "field": "name", "cursor_pos": 3}
    )
    row1.name_field.update = MagicMock()
    row1.name_field.focus = MagicMock()
    row1.seq_field = ft.TextField(
        value="GGCC", data={"idx": 1, "field": "seq", "cursor_pos": 0}
    )
    row1.seq_field.update = MagicMock()
    row1.seq_field.focus = MagicMock()

    ctrl.input_view.primer_input.primers_list.controls = [row0, row1]

    # Tab navigation
    ctrl.input_view._currently_focused_control = row0.name_field
    ev_tab = MagicMock()
    ev_tab.key = "Tab"
    handle_keyboard_event(ctrl, ev_tab)

    ctrl.input_view._currently_focused_control = row0.seq_field
    handle_keyboard_event(ctrl, ev_tab)

    ctrl.input_view._currently_focused_control = row1.seq_field
    handle_keyboard_event(ctrl, ev_tab)

    # Arrow keys navigation
    ctrl.input_view._currently_focused_control = row0.name_field
    ev_right = MagicMock()
    ev_right.key = "Arrow Right"
    handle_keyboard_event(ctrl, ev_right)

    ctrl.input_view._currently_focused_control = row0.seq_field
    ev_left = MagicMock()
    ev_left.key = "Arrow Left"
    handle_keyboard_event(ctrl, ev_left)

    ev_down = MagicMock()
    ev_down.key = "Arrow Down"
    handle_keyboard_event(ctrl, ev_down)

    ev_up = MagicMock()
    ev_up.key = "Arrow Up"
    handle_keyboard_event(ctrl, ev_up)

    # 6. System SHA & lifecycle helpers
    with (
        patch("subprocess.run", side_effect=subprocess.SubprocessError("fail")),
        patch("os.path.exists", side_effect=lambda p: ".git-sha" in p),
        patch("builtins.open", mock_open(read_data="dist_sha_12345")),
    ):
        assert _get_sha(full=False) == "dist_sha_12345"
        assert _get_sha(full=True) == "dist_sha_12345"

    with (
        patch("subprocess.run", side_effect=subprocess.SubprocessError("fail")),
        patch("os.path.exists", return_value=False),
    ):
        assert _get_sha(full=False) == "unknown"

    with patch.dict("sys.modules", {"amplifyp": None}):
        assert get_version() is not None

    ctrl._confirm_dialog = MagicMock(open=True)  # type: ignore[assignment]
    confirm_dismiss(ctrl, MagicMock())
    assert cast(Any, ctrl._confirm_dialog).open is False

    mock_page.window.destroy = AsyncMock()

    await confirm_exit_async(ctrl)
    confirm_exit(ctrl, MagicMock())

    # Screenshots & delayed quit
    mock_page.take_screenshot = AsyncMock(return_value=b"png_bytes")
    scr_path = await capture_view_screenshot_async(
        mock_page, "test_view", screenshots_dir=tmp_path / "test_screens"
    )
    assert scr_path.exists()

    ctrl.auto_close = True
    ctrl.export_screenshots = True
    await capture_all_views_async(ctrl)

    ctrl.state_file = "state.yaml"
    with patch.object(ctrl, "_restore_state_from_file"):
        await restore_state_and_auto_close_async(ctrl)

    await auto_close_and_quit_delayed(ctrl)

    ev_win = MagicMock()
    ev_win.data = "close"
    on_window_event(ctrl, ev_win)
    assert ctrl._confirm_dialog is not None


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_data_helpers_and_system_utilities_extra() -> None:
    """Test edge branches in data helpers and system lifecycle."""
    import asyncio
    from importlib.metadata import PackageNotFoundError
    from pathlib import Path
    from typing import Any
    from unittest.mock import AsyncMock, mock_open, patch

    from amplifyp.gui.controller import GUIController
    from amplifyp.gui.utils.data_helpers import (
        apply_parsed_state,
        confirm_clear,
        dismiss_clear,
        load_last_state,
        load_state,
        restore_state_from_file,
        save_state,
    )
    from amplifyp.gui.utils.gui_helpers import (
        focus_async,
        handle_keyboard_event,
    )
    from amplifyp.gui.utils.system import (
        _get_sha,
        auto_close_and_quit_delayed,
        capture_all_views_async,
        capture_view_screenshot_async,
        confirm_exit_async,
        get_version,
        on_window_event,
        restore_state_and_auto_close_async,
    )
    from amplifyp.gui.views.input.primer.row import PrimerRow

    mock_page = MagicMock(spec=ft.Page)
    mock_page.services = []
    mock_page.overlay = []
    mock_page.web = False
    mock_page.window = MagicMock()
    mock_page.window.destroy = AsyncMock(
        side_effect=RuntimeError("Already destroyed")
    )
    mock_page.controls = []

    ctrl = GUIController(page=mock_page)
    ctrl.initialise()

    # 1. focus_async
    async def dummy_coro() -> int:
        return 42

    await focus_async(dummy_coro())

    # 2. data_helpers missing lines
    # load_last_state valid desktop yaml
    with (
        patch.object(Path, "exists", return_value=True),
        patch(
            "builtins.open",
            mock_open(read_data="template: ATGC\nprimers: []\n"),
        ),
    ):
        load_last_state(ctrl)

    # restore_state_from_file valid
    with patch(
        "builtins.open",
        mock_open(
            read_data="input: {template: ATGC}\nsettings: {dark_mode: true}\n"
        ),
    ):
        restore_state_from_file(ctrl, "valid.yaml")

    # apply_parsed_state without input and with settings
    apply_parsed_state(
        ctrl,
        {"template": "ATGC", "settings": {"dark_mode": False}},
        ignore_settings=False,
    )

    # save_state error
    with patch(
        "amplifyp.gui.utils.data_helpers.save_and_write_file",
        side_effect=ValueError("save fail"),
    ):
        await save_state(ctrl, MagicMock())

    # load_state pick returns None and yaml error
    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        return_value=None,
    ):
        await load_state(ctrl, MagicMock())
    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        return_value=": [invalid yaml",
    ):
        await load_state(ctrl, MagicMock())

    # confirm_clear & dismiss_clear with _clear_dialog
    ctrl._clear_dialog = MagicMock()  # type: ignore[assignment]
    confirm_clear(ctrl, MagicMock())
    dismiss_clear(ctrl, MagicMock())

    # 3. gui_helpers keyboard missing lines
    # Missing idx or field data (line 232)
    ctrl.input_view._currently_focused_control = ft.TextField(
        data={"missing_idx": 1}
    )
    handle_keyboard_event(ctrl, MagicMock(key="Tab"))

    # Non-textfield focused (line 235)
    ctrl.input_view._currently_focused_control = ft.Container(
        data={"idx": 0, "field": "name"}
    )
    handle_keyboard_event(ctrl, MagicMock(key="Tab"))

    # Tab on name with empty row list (line 256)
    ctrl.input_view.primer_input.primers_list.controls = []
    ctrl.input_view._currently_focused_control = ft.TextField(
        data={"idx": 0, "field": "name"}
    )
    handle_keyboard_event(ctrl, MagicMock(key="Tab"))

    # Tab on unknown field (line 272)
    ctrl.input_view._currently_focused_control = ft.TextField(
        data={"idx": 0, "field": "unknown_field"}
    )
    handle_keyboard_event(ctrl, MagicMock(key="Tab"))

    # Arrow Right with cursor not at end (line 284)
    ctrl.input_view._currently_focused_control = ft.TextField(
        value="Fwd", data={"idx": 0, "field": "name", "cursor_pos": 1}
    )
    handle_keyboard_event(ctrl, MagicMock(key="Arrow Right"))

    # Arrow Right with empty row list (line 291)
    ctrl.input_view._currently_focused_control = ft.TextField(
        value="Fwd", data={"idx": 0, "field": "name", "cursor_pos": 3}
    )
    handle_keyboard_event(ctrl, MagicMock(key="Arrow Right"))

    # Arrow Left with cursor_pos != 0 (line 299)
    ctrl.input_view._currently_focused_control = ft.TextField(
        value="ATGC", data={"idx": 0, "field": "seq", "cursor_pos": 2}
    )
    handle_keyboard_event(ctrl, MagicMock(key="Arrow Left"))

    # Arrow Left on seq with empty row list (line 306)
    ctrl.input_view._currently_focused_control = ft.TextField(
        value="ATGC", data={"idx": 0, "field": "seq", "cursor_pos": 0}
    )
    handle_keyboard_event(ctrl, MagicMock(key="Arrow Left"))

    # Arrow Left on name field (line 313)
    ctrl.input_view._currently_focused_control = ft.TextField(
        value="Fwd", data={"idx": 0, "field": "name", "cursor_pos": 0}
    )
    handle_keyboard_event(ctrl, MagicMock(key="Arrow Left"))

    # Key other than Arrow Up/Down (line 322)
    handle_keyboard_event(ctrl, MagicMock(key="Enter"))

    tasks = []

    def mock_run_task(fn: Any, *args: Any) -> None:
        task = asyncio.create_task(fn(*args))
        tasks.append(task)

    mock_page.run_task = mock_run_task

    row0 = MagicMock(spec=PrimerRow)
    row0.idx = 0
    row0.checkbox = MagicMock(value=True)

    async def async_focus() -> None:
        pass

    row0.name_field = ft.TextField(
        value="Fwd", data={"idx": 0, "field": "name", "cursor_pos": 3}
    )
    row0.name_field.update = MagicMock()
    row0.name_field.focus = MagicMock(side_effect=lambda: async_focus())
    row0.seq_field = ft.TextField(
        value="ATGC", data={"idx": 0, "field": "seq", "cursor_pos": 0}
    )
    row0.seq_field.update = MagicMock()
    row0.seq_field.focus = MagicMock(side_effect=lambda: async_focus())

    row1 = MagicMock(spec=PrimerRow)
    row1.idx = 1
    row1.checkbox = MagicMock(value=True)
    row1.name_field = ft.TextField(
        value="Rev", data={"idx": 1, "field": "name", "cursor_pos": 3}
    )
    row1.name_field.update = MagicMock(
        side_effect=RuntimeError("simulated error")
    )
    row1.name_field.focus = MagicMock(side_effect=lambda: async_focus())
    row1.seq_field = ft.TextField(
        value="GGCC", data={"idx": 1, "field": "seq", "cursor_pos": 0}
    )
    row1.seq_field.update = MagicMock()
    row1.seq_field.focus = MagicMock(side_effect=lambda: async_focus())

    # Restore rows
    ctrl.input_view.primer_input.primers_list.controls = [row0, row1]

    # Tab with focus coroutine
    ctrl.input_view._currently_focused_control = row0.name_field
    ev_tab = MagicMock(key="Tab")
    handle_keyboard_event(ctrl, ev_tab)

    # Arrow Right with focus coroutine
    handle_keyboard_event(ctrl, MagicMock(key="Arrow Right"))

    # Arrow Left on seq field with focus coroutine
    ctrl.input_view._currently_focused_control = row0.seq_field
    handle_keyboard_event(ctrl, MagicMock(key="Arrow Left"))

    # Arrow Down with update RuntimeError and focus coroutine
    ctrl.input_view._currently_focused_control = row0.name_field
    handle_keyboard_event(ctrl, MagicMock(key="Arrow Down"))

    if tasks:
        await asyncio.gather(*tasks)

    # 4. system.py missing lines
    # OSError in git HEAD read
    with (
        patch("subprocess.run", side_effect=subprocess.SubprocessError("fail")),
        patch("os.path.exists", return_value=True),
        patch("builtins.open", side_effect=OSError("disk error")),
    ):
        assert _get_sha(full=False) == "unknown"

    # OSError in .git-sha read
    with (
        patch("subprocess.run", side_effect=subprocess.SubprocessError("fail")),
        patch("os.path.exists", side_effect=lambda p: ".git-sha" in p),
        patch("builtins.open", side_effect=OSError("read fail")),
    ):
        assert _get_sha(full=False) == "unknown"

    # PackageNotFoundError in get_version
    with (
        patch.dict("sys.modules", {"amplifyp": None}),
        patch("importlib.metadata.version", side_effect=PackageNotFoundError),
    ):
        assert "unknown" in get_version()

    await confirm_exit_async(ctrl)

    # Window event with WindowEventType.CLOSE and auto_reload_on_startup False
    ctrl.settings["auto_reload_on_startup"] = False
    ctrl._confirm_dialog = None
    ev_close = MagicMock(data="other", type=ft.WindowEventType.CLOSE)
    on_window_event(ctrl, ev_close)
    assert ctrl._confirm_dialog is not None
    assert "Unsaved changes" in ctrl._confirm_dialog.content.value
    # Second invocation with existing dialog
    on_window_event(ctrl, ev_close)

    # capture_view_screenshot_async failure
    mock_page.take_screenshot = AsyncMock(return_value=None)
    with pytest.raises(RuntimeError):
        await capture_view_screenshot_async(mock_page, "empty_view")

    # capture_all_views_async with enabled buttons and auto_close
    mock_page.take_screenshot = AsyncMock(return_value=b"bytes")
    ctrl.pcr_view.open_all_cards = MagicMock()
    ctrl.pcr_view.run_pcr = MagicMock()
    ctrl.dimers_view.run_analysis = MagicMock()

    with patch.object(ctrl, "update_pcr_button_state"):
        pcr_btn = MagicMock(disabled=False)
        ctrl.pcr_button_ref.current = pcr_btn
        dimers_btn = MagicMock(disabled=False)
        ctrl.dimers_button_ref.current = dimers_btn
        ctrl.auto_close = True
        await capture_all_views_async(ctrl)
        ctrl.pcr_view.open_all_cards.assert_called_once()
        ctrl.pcr_view.run_pcr.assert_called_once()
        ctrl.dimers_view.run_analysis.assert_called_once()

    # capture_all_views_async exception with confirm_exit_async
    # raising RuntimeError
    with (
        patch(
            "amplifyp.gui.utils.system.capture_view_screenshot_async",
            side_effect=Exception("scr fail"),
        ),
        patch.object(
            ctrl, "confirm_exit_async", side_effect=RuntimeError("closed")
        ),
    ):
        await capture_all_views_async(ctrl)

    # auto_close_and_quit_delayed exception with confirm_exit_async
    # raising RuntimeError
    with (
        patch.object(
            ctrl.pcr_view, "run_pcr", side_effect=Exception("pcr fail")
        ),
        patch.object(
            ctrl, "confirm_exit_async", side_effect=RuntimeError("closed")
        ),
    ):
        await auto_close_and_quit_delayed(ctrl)

    # restore_state_and_auto_close_async with auto_close and state_file
    ctrl.export_screenshots = False
    ctrl.auto_close = True
    ctrl.state_file = "state.yaml"
    with (
        patch.object(ctrl, "_restore_state_from_file"),
        patch.object(
            ctrl, "_auto_close_and_quit_delayed", new_callable=AsyncMock
        ),
    ):
        await restore_state_and_auto_close_async(ctrl)
