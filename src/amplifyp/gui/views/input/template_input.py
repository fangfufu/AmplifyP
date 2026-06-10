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


class TemplateInput(ft.Container):  # type: ignore[misc]
    """Input component for DNA template sequence."""

    def __init__(
        self,
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
        self.settings = settings
        self.input_data = input_data

        font_family = self.settings.get("font_family", "Roboto Mono")
        self.template_sequence = ft.TextField(
            dense=True,
            multiline=True,
            expand=True,
            hint_text="Enter DNA sequence here...",
            content_padding=0,
            on_change=on_change_handler,
            on_focus=handle_field_focus,
            on_blur=handle_field_blur,
            on_submit=handle_field_submit,
            text_style=ft.TextStyle(font_family=font_family),
        )
        self.template_circular = ft.Checkbox(
            label="Circular Template",
            value=False,
            on_change=on_change_handler,
        )
        self.circular_container = ft.Container(
            content=self.template_circular,
            border=ft.Border.all(1, GUIColors.OUTLINE),
            border_radius=5,
            padding=ft.Padding(10, 0, 10, 0),
            height=32,
            alignment=ft.Alignment(0, 0),
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
                ft.Row(
                    [
                        ft.Text("Template Sequence", weight=ft.FontWeight.BOLD),
                        ft.Row(
                            [
                                self.circular_container,
                                self.clear_template_button,
                            ],
                            spacing=10,
                        ),
                    ],
                    alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                    height=40,
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
