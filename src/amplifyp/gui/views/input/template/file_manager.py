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

"""File management helpers for DNA template sequence input."""

from __future__ import annotations

from typing import TYPE_CHECKING

import flet as ft

if TYPE_CHECKING:
    from .input import TemplateInput


async def load_template_click(
    template_input: TemplateInput, e: ft.Event
) -> None:
    """Open file picker to load template sequence from a TXT file."""
    from amplifyp.gui.utils.data_helpers import pick_and_read_file

    content = await pick_and_read_file(
        page=template_input.app_page,
        dialog_title="Load",
        allowed_extensions=["txt"],
        show_notification=template_input._show_notification,
    )
    if content is None:
        return

    template_input.input_data.template = content
    template_input.update_ui()
    template_input.on_change_handler(None)
    template_input._show_notification("Template loaded successfully.")


async def save_template_click(
    template_input: TemplateInput, e: ft.Event
) -> None:
    """Save template sequence to a TXT file."""
    template_content = template_input.input_data.template
    if not template_content.strip():
        template_input._show_notification("No template to save.")
        return

    from amplifyp.gui.utils.data_helpers import save_and_write_file

    await save_and_write_file(
        page=template_input.app_page,
        dialog_title="Save",
        file_name="template.txt",
        allowed_extensions=["txt"],
        content=template_content,
        show_notification=template_input._show_notification,
        success_message_desktop="Template saved successfully.",
        success_message_web="Template ready for download!",
    )
