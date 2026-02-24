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

"""Tests for GUI state saving and loading in Web (Pyodide) environment."""

import sys
from collections.abc import Generator
from unittest.mock import MagicMock, patch

import flet as ft
import pytest


@pytest.fixture  # type: ignore[untyped-decorator]
def mock_web_modules() -> Generator[tuple[MagicMock, MagicMock], None, None]:
    """Mock js and pyodide modules."""
    mock_js = MagicMock()
    mock_pyodide = MagicMock()
    mock_pyodide.ffi.to_js = lambda x, **kwargs: x
    mock_pyodide.ffi.create_proxy = lambda x: x

    with patch.dict(
        sys.modules,
        {
            "js": mock_js,
            "pyodide": mock_pyodide,
            "pyodide.ffi": mock_pyodide.ffi,
        },
    ):
        yield mock_js, mock_pyodide


@pytest.mark.anyio  # type: ignore[untyped-decorator]
async def test_web_state_save_load(
    mock_web_modules: tuple[MagicMock, MagicMock],
) -> None:
    """Test saving and loading state in a mock Pyodide environment."""
    mock_js, _ = mock_web_modules

    # Import main inside the patched environment
    from amplifyp.gui.views import InputView
    from main import main

    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = True

    # We need to capture the views created in main
    # main() creates views and calls page.add(container)
    # The container content is initially input_view.

    # Run main logic
    with patch("main._is_pyodide", return_value=True):
        main(mock_page)

        # Extract the container added to the page
        assert mock_page.add.called
        container = mock_page.add.call_args[0][0]
        assert isinstance(container, ft.Container)

        # Extract views.
        # main.py defines:
        # input_view = InputView(...)
        # settings_view = SettingsView(page)
        # result_view = ResultView(page, input_view, settings_view)
        # It doesn't expose them directly, but we can access them through the
        # closures of the event handlers or by inspecting the container content.

        # Initially, view_container.content = input_view
        input_view = container.content
        assert isinstance(input_view, InputView)

        # Get actions from appbar to find buttons
        actions = mock_page.appbar.actions
        assert actions

        # Helper to find button by icon
        def find_button_by_icon(icon_name: str) -> ft.IconButton | None:
            for control in actions:
                if (
                    isinstance(control, ft.IconButton)
                    and control.icon == icon_name
                ):
                    return control
            return None

        save_btn = find_button_by_icon(ft.Icons.SAVE)
        load_btn = find_button_by_icon(ft.Icons.UPLOAD_FILE)

        assert save_btn
        assert load_btn

        # --- TEST SAVE ---

        # Set some state
        input_view.template_sequence.value = "ACGT"
        input_view.template_circular.value = True

        # To get the settings_view instance, we can look at save_state closure
        # if we could, but simpler is to simulate the save and see what it
        # grabs.
        # main.py: state["settings"] = settings_view.get_state()

        # Trigger Save
        await save_btn.on_click(MagicMock(spec=ft.ControlEvent))

        # Verify js.postMessage called
        assert mock_js.postMessage.called
        call_args = mock_js.postMessage.call_args[0][0]

        assert call_args["type"] == "save_file"
        assert call_args["filename"] == "amplify_gui_state.yaml"
        assert "ACGT" in call_args["content"]
        assert "template_circular: true" in call_args["content"]

        # --- TEST LOAD ---

        # Reset Mock
        mock_js.postMessage.reset_mock()

        # Trigger Load
        await load_btn.on_click(MagicMock(spec=ft.ControlEvent))

        # Verify it asked to open file picker
        mock_js.postMessage.assert_called_with("open_file_picker")

        # Verify callback was registered
        # main.py: js.self.custom_file_callback = create_proxy(on_file_content)
        # Note: js.self is mock_js.self
        assert mock_js.self.custom_file_callback
        callback = mock_js.self.custom_file_callback

        # Create mocked YAML content to load
        yaml_content = """
template: GGGCCC
template_circular: false
primers: []
settings:
  primability_cutoff: 0.5
"""
        # Execute callback
        callback(yaml_content)

        # Verify state updated
        # We need to access the input_view and settings_view to verify.
        # input_view we already have.
        assert input_view.template_sequence.value == "GGGCCC"
        assert not input_view.template_circular.value

        # To verify settings_view, we can try to save again and check the
        # output.
        mock_js.postMessage.reset_mock()
        await save_btn.on_click(MagicMock(spec=ft.ControlEvent))

        save_args = mock_js.postMessage.call_args[0][0]
        assert (
            "primability_cutoff: '0.5'" in save_args["content"]
            or "primability_cutoff: 0.5" in save_args["content"]
        )
