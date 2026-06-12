# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""Input component for DNA template sequence."""

from typing import Any

import flet as ft

from amplifyp.gui.settings import GUIColors, GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import clean_sequence

from .template_file_manager import TemplateFileManager


class TemplateInput(ft.Container):  # type: ignore[misc]
    """Input component for DNA template sequence."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        input_data: GUIInput,
        on_change_handler: Any,
        handle_field_focus: Any,
        handle_field_blur: Any,
        handle_field_submit: Any,
        clear_template_callback: Any,
    ) -> None:
        """Initialize the TemplateInput component."""
        super().__init__(expand=5)
        self.app_page = page
        self.settings = settings
        self.input_data = input_data

        # Template File Manager
        self.file_manager = TemplateFileManager(
            page=self.app_page,
            input_data=self.input_data,
            on_update_ui=self.update_ui,
            on_change_handler=on_change_handler,
            show_snackbar=self._show_snackbar,
        )

        font_family = self.settings.get("font_family", "Roboto Mono")
        self.template_sequence = ft.TextField(
            dense=True,
            multiline=True,
            expand=True,
            hint_text="Enter DNA sequence here...",
            border=ft.InputBorder.NONE,
            content_padding=10,
            on_change=on_change_handler,
            on_focus=handle_field_focus,
            on_blur=handle_field_blur,
            on_submit=handle_field_submit,
            text_style=ft.TextStyle(font_family=font_family),
        )
        self.template_circular = ft.Checkbox(
            label="Circular",
            value=False,
            on_change=on_change_handler,
        )
        self.circular_container = ft.Container(
            content=self.template_circular,
            border=ft.Border.all(1, GUIColors.OUTLINE),
            border_radius=5,
            padding=ft.Padding(0, 0, 10, 0),
            height=32,
            width=110,
        )

        self.save_template_button = ft.FilledTonalButton(
            "Save",
            icon=ft.Icons.FILE_DOWNLOAD,
            on_click=self._save_template_click,
            tooltip="Save template to TXT",
            height=32,
        )

        self.load_template_button = ft.FilledTonalButton(
            "Load",
            icon=ft.Icons.FILE_OPEN,
            on_click=self._load_template_click,
            tooltip="Load template from TXT",
            height=32,
        )

        self.clear_template_button = ft.OutlinedButton(
            "Clear",
            icon=ft.Icons.DELETE_OUTLINE,
            tooltip="Clear Template",
            on_click=clear_template_callback,
            height=32,
        )

        self.content = ft.Column(
            [
                ft.ResponsiveRow(
                    [
                        ft.Container(
                            content=ft.Row(
                                [
                                    self.circular_container,
                                    ft.Text(
                                        "Template Sequence",
                                        weight=ft.FontWeight.BOLD,
                                        no_wrap=True,
                                    ),
                                ],
                                spacing=10,
                                tight=True,
                                wrap=True,
                            ),
                            col={
                                "xl": 5,
                                "lg": 12,
                                "md": 12,
                                "sm": 12,
                                "xs": 12,
                            },
                        ),
                        ft.Container(
                            content=ft.Row(
                                [
                                    self.load_template_button,
                                    self.save_template_button,
                                    self.clear_template_button,
                                ],
                                spacing=10,
                                tight=True,
                                wrap=True,
                            ),
                            col={
                                "xl": 7,
                                "lg": 12,
                                "md": 12,
                                "sm": 12,
                                "xs": 12,
                            },
                            alignment=ft.Alignment(1, 0),
                        ),
                    ],
                ),
                ft.Container(
                    content=ft.ListView(
                        [self.template_sequence],
                        expand=True,
                        scroll=ft.ScrollMode.ALWAYS,
                    ),
                    expand=True,
                    border=ft.Border.all(1, GUIColors.OUTLINE),
                    border_radius=5,
                    padding=0,
                ),
            ],
            expand=True,
            spacing=5,
        )

    async def _load_template_click(self, e: ft.ControlEvent) -> None:
        """Open file picker to load template sequence."""
        await self.file_manager.load_template_click(e)

    async def _save_template_click(self, e: ft.ControlEvent) -> None:
        """Save template sequence."""
        await self.file_manager.save_template_click(e)

    def _show_snackbar(self, message: str) -> None:
        """Show a snackbar message."""
        if not hasattr(self, "_snack_bar"):
            self._snack_bar = ft.SnackBar(ft.Text(""), open=False)
            self.app_page.overlay.append(self._snack_bar)
        self._snack_bar.content = ft.Text(message)
        self._snack_bar.open = True
        self.app_page.update()

    def sync_to_state(self) -> None:
        """Sync template text field to the central state."""
        self.template_sequence.value = self.template_sequence.value or ""
        self.input_data.template = clean_sequence(
            str(self.template_sequence.value)
        )
        self.input_data.template_circular = bool(self.template_circular.value)

    def update_ui(self) -> None:
        """Update template UI elements to match central state."""
        font_family = self.settings.get("font_family", "Roboto Mono")
        self.template_sequence.text_style = ft.TextStyle(
            font_family=font_family
        )
        self.template_sequence.value = self.input_data.template
        self.template_circular.value = self.input_data.template_circular
