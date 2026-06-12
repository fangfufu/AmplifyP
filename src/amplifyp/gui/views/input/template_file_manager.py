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
        try:
            files = await ft.FilePicker().pick_files(
                dialog_title="Load",
                allowed_extensions=["txt"],
                file_type=ft.FilePickerFileType.CUSTOM,
                with_data=True,
            )
            if not files:
                return

            file = files[0]
            if file.bytes is not None:
                content = file.bytes.decode("utf-8")
            else:
                if not file.path:
                    self.show_snackbar("Error: Could not read file content.")
                    return
                with open(file.path, encoding="utf-8") as f:
                    content = f.read()

            self.input_data.template = content
            self.on_update_ui()
            if self.on_change_handler:
                self.on_change_handler(None)
            self.show_snackbar("Template loaded successfully.")

        except Exception as ex:
            self.show_snackbar(f"Error loading file: {ex}")

    async def save_template_click(self, e: ft.ControlEvent) -> None:
        """Save template sequence to a TXT file."""
        template_content = self.input_data.template
        if not template_content.strip():
            self.show_snackbar("No template to save.")
            return

        try:
            file_path = await ft.FilePicker().save_file(
                dialog_title="Save",
                file_name="template.txt",
                allowed_extensions=["txt"],
                file_type=ft.FilePickerFileType.CUSTOM,
                src_bytes=template_content.encode("utf-8"),
            )
            if self.app_page.web:
                self.show_snackbar("Template ready for download!")
            else:
                if file_path is None:
                    return
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(template_content)
                self.show_snackbar("Template saved successfully.")
        except Exception as ex:
            self.show_snackbar(f"Error saving file: {ex}")
