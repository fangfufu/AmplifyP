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

"""Tests for GUI state saving and loading in Web environment."""

from typing import Any
from unittest.mock import AsyncMock, MagicMock, patch

import flet as ft
import pytest


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_web_state_save_load() -> None:
    """Test saving and loading state using native FilePicker on Web."""
    # Import main inside the test
    from amplifyp.gui.views import InputView
    from main import main

    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = True

    # Mock FilePicker class
    mock_file_picker_instance = MagicMock(spec=ft.FilePicker)
    mock_file_picker_instance.save_file = AsyncMock(
        return_value="amplify_gui_state.yaml"
    )

    with patch(
        "amplifyp.gui.app.ft.FilePicker", return_value=mock_file_picker_instance
    ):
        # Run main logic
        main(mock_page)

        # Extract the container added to the page
        assert mock_page.add.called
        container = mock_page.add.call_args[0][0]
        assert isinstance(container, ft.Container)

        # Extract input_view
        input_view = container.content
        assert isinstance(input_view, InputView)

        # Get actions from appbar to find buttons
        actions = mock_page.appbar.actions
        assert actions

        # Helper to find button by icon
        def find_button_by_icon(icon_name: str) -> Any:
            for control in actions:
                if hasattr(control, "icon") and control.icon == icon_name:
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

        # Trigger Save
        await save_btn.on_click(MagicMock(spec=ft.ControlEvent))

        # Verify save_file called with the correct parameters
        mock_file_picker_instance.save_file.assert_called_once()
        kwargs = mock_file_picker_instance.save_file.call_args[1]
        assert kwargs["file_name"] == "amplify_gui_state.yaml"
        assert kwargs["file_type"] == ft.FilePickerFileType.CUSTOM
        assert kwargs["allowed_extensions"] == ["yaml", "yml"]

        # Verify src_bytes content
        src_bytes = kwargs["src_bytes"]
        assert b"ACGT" in src_bytes
        assert b"template_circular: true" in src_bytes

        # --- TEST LOAD ---
        # Create a mock FilePickerFile with yaml content
        yaml_content = """
template: GGGCCC
template_circular: false
primers: []
settings:
  primability_cutoff: 0.5
"""
        mock_file = MagicMock(spec=ft.FilePickerFile)
        mock_file.name = "amplify_gui_state.yaml"
        mock_file.bytes = yaml_content.encode("utf-8")
        mock_file.path = None

        # Mock pick_files to return the mock file
        mock_file_picker_instance.pick_files = AsyncMock(
            return_value=[mock_file]
        )

        # Trigger Load
        await load_btn.on_click(MagicMock(spec=ft.ControlEvent))

        # Verify pick_files was called
        mock_file_picker_instance.pick_files.assert_called_once_with(
            dialog_title="Load State",
            allowed_extensions=["yaml", "yml"],
            file_type=ft.FilePickerFileType.CUSTOM,
            with_data=True,
        )

        # Verify state updated in input_view
        assert input_view.template_sequence.value == "GGGCCC"
        assert not input_view.template_circular.value
