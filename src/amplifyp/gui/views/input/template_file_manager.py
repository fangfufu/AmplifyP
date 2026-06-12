# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Handles template sequence file loading and saving."""

from typing import Any

import flet as ft

from amplifyp.gui.user_data import GUIInput


class TemplateFileManager:
    """Handles template sequence file operations."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput,
        on_update_ui: Any,
        on_change_handler: Any,
        show_snackbar: Any,
    ) -> None:
        """Initialize the TemplateFileManager."""
        self.app_page = page
        self.input_data = input_data
        self.on_update_ui = on_update_ui
        self.on_change_handler = on_change_handler
        self.show_snackbar = show_snackbar

    async def load_template_click(self, e: ft.ControlEvent) -> None:
        """Open file picker to load template sequence from a TXT file."""
        from amplifyp.gui.util import pick_and_read_file

        content = await pick_and_read_file(
            page=self.app_page,
            dialog_title="Load",
            allowed_extensions=["txt"],
            show_snackbar=self.show_snackbar,
        )
        if content is None:
            return

        self.input_data.template = content
        self.on_update_ui()
        if self.on_change_handler:
            self.on_change_handler(None)
        self.show_snackbar("Template loaded successfully.")

    async def save_template_click(self, e: ft.ControlEvent) -> None:
        """Save template sequence to a TXT file."""
        template_content = self.input_data.template
        if not template_content.strip():
            self.show_snackbar("No template to save.")
            return

        from amplifyp.gui.util import save_and_write_file

        await save_and_write_file(
            page=self.app_page,
            dialog_title="Save",
            file_name="template.txt",
            allowed_extensions=["txt"],
            content=template_content,
            show_snackbar=self.show_snackbar,
            success_message_desktop="Template saved successfully.",
            success_message_web="Template ready for download!",
        )
