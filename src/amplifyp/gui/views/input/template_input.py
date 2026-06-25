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

"""Input component for DNA template sequence."""

from __future__ import annotations

from collections.abc import Callable

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import NotificationHelper, clean_sequence


class TemplateInput(ft.Container):  # type: ignore[misc]
    """Input component for DNA template sequence."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        input_data: GUIInput,
        on_change_handler: Callable[[ft.Event | None], None],
        handle_field_focus: Callable[[ft.Event[ft.TextField]], None],
        handle_field_blur: Callable[[ft.Event[ft.TextField]], None],
        handle_field_submit: Callable[[ft.Event[ft.TextField]], None],
        clear_template_callback: Callable[[ft.Event | None], None],
    ) -> None:
        """Initialise the TemplateInput component."""
        super().__init__(expand=5)
        self.app_page = page
        self.settings = settings
        self.input_data = input_data
        self.on_change_handler = on_change_handler

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
            border=ft.Border.all(1, GUIColours.OUTLINE),
            border_radius=5,
            padding=ft.Padding(0, 0, 0, 0),
            height=32,
            width=100,
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
                    border=ft.Border.all(1, GUIColours.OUTLINE),
                    border_radius=5,
                    padding=0,
                ),
            ],
            expand=True,
            spacing=5,
        )

    async def _load_template_click(self, e: ft.Event) -> None:
        """Open file picker to load template sequence from a TXT file.

        Args:
            e: The Flet control event triggered by the load button click.
        """
        from amplifyp.gui.util import pick_and_read_file

        content = await pick_and_read_file(
            dialog_title="Load",
            allowed_extensions=["txt"],
            show_notification=self._show_notification,
        )
        if content is None:
            return

        self.input_data.template = content
        self.update_ui()
        self.on_change_handler(None)
        self._show_notification("Template loaded successfully.")

    async def _save_template_click(self, e: ft.Event) -> None:
        """Save template sequence to a TXT file.

        Args:
            e: The Flet control event triggered by the save button click.
        """
        template_content = self.input_data.template
        if not template_content.strip():
            self._show_notification("No template to save.")
            return

        from amplifyp.gui.util import save_and_write_file

        await save_and_write_file(
            page=self.app_page,
            dialog_title="Save",
            file_name="template.txt",
            allowed_extensions=["txt"],
            content=template_content,
            show_notification=self._show_notification,
            success_message_desktop="Template saved successfully.",
            success_message_web="Template ready for download!",
        )

    def _show_notification(self, message: str) -> None:
        """Show a notification message.

        Args:
            message: The message to display in the notification.
        """
        if not hasattr(self, "_notification_helper"):
            self._notification_helper = NotificationHelper(self.app_page)
        self._notification_helper.show_message(message)

    def sync_to_state(self) -> None:
        """Sync template text field to the central state.

        Reads the current UI values and writes them into the central
        ``GUIInput`` state object.
        """
        self.template_sequence.value = self.template_sequence.value or ""
        self.input_data.template = clean_sequence(
            str(self.template_sequence.value)
        )
        self.input_data.template_circular = bool(self.template_circular.value)

    def update_ui(self) -> None:
        """Update template UI elements to match central state.

        Applies values from the central ``GUIInput`` state object to the
        template text field and circular checkbox controls.
        """
        font_family = self.settings.get("font_family", "Roboto Mono")
        self.template_sequence.text_style = ft.TextStyle(
            font_family=font_family
        )
        self.template_sequence.value = self.input_data.template
        self.template_circular.value = self.input_data.template_circular
