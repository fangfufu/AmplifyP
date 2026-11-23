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

"""Shared base class for settings weight table components."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours


class WeightTableBase(ft.Column):  # type: ignore[misc]
    """Base Flet Column for labelled, bordered DataTable weight editors."""

    def __init__(
        self,
        label: str,
        prefix: str,
        settings_map: dict[str, Any],
        font_size_default: int,
        font_size_table_header: int,
        on_change_handler: Callable[[ft.Event | None], None] | None = None,
    ) -> None:
        """Store shared configuration for weight table subclasses."""
        self.label = label
        self.prefix = prefix
        self.settings_map = settings_map
        self.font_size_default = font_size_default
        self.font_size_table_header = font_size_table_header
        self.on_change_handler = on_change_handler

    def _wrap_centred(self, control: ft.Control, width: int) -> ft.Container:
        """Wrap a control in a fixed-width, centred Container."""
        return ft.Container(
            content=control,
            width=width,
            alignment=ft.Alignment(0, 0),
        )

    def _build_header_label(self, text: str, width: int) -> ft.Container:
        """Create a centred, bold table header or row-label Container."""
        return self._wrap_centred(
            ft.Text(
                text,
                weight=ft.FontWeight.BOLD,
                size=self.font_size_table_header,
            ),
            width=width,
        )

    def _style_field(self, field: ft.TextField, width: int) -> ft.TextField:
        """Apply borderless, centred styling to a table cell TextField."""
        field.border = ft.NoInputBorder()
        field.width = width
        field.height = 30
        field.content_padding = 0
        field.dense = True
        field.text_align = ft.TextAlign.CENTER
        field.text_vertical_align = ft.VerticalAlignment.CENTER
        return field

    def _create_data_table(
        self,
        columns: list[ft.DataColumn],
        rows: list[ft.DataRow],
        *,
        column_spacing: int,
        horizontal_margin: int,
    ) -> ft.DataTable:
        """Create a styled DataTable matching settings view conventions."""
        return ft.DataTable(
            border=ft.Border.all(1, GUIColours.TRANSPARENT),
            vertical_lines=ft.BorderSide(1, GUIColours.DIVIDER_GREY),
            horizontal_lines=ft.BorderSide(1, GUIColours.DIVIDER_GREY),
            column_spacing=column_spacing,
            horizontal_margin=horizontal_margin,
            heading_row_color=GUIColours.INFO_HEADER_BG,
            heading_row_height=40,
            data_row_min_height=36,
            data_row_max_height=36,
            columns=columns,
            rows=rows,
        )

    def _create_table_container(
        self, table: ft.DataTable, width: int | None = None
    ) -> ft.Container:
        """Wrap a DataTable in a rounded outline Container."""
        return ft.Container(
            content=table,
            border=ft.Border.all(1, GUIColours.OUTLINE),
            border_radius=5,
            padding=0,
            width=width,
        )

    def _create_reset_button(
        self, on_click: Callable[[ft.Event | None], None]
    ) -> ft.TextButton:
        """Create a standard Reset button for the table."""
        return ft.TextButton(
            "Reset",
            icon=ft.Icons.RESTORE,
            on_click=on_click,
        )

    def _init_column(
        self,
        table_control: ft.Control,
        action_controls: list[ft.Control] | None = None,
        *,
        spacing: int,
    ) -> None:
        """Initialise the Flet Column with title, table, and actions."""
        controls: list[ft.Control] = [
            ft.Text(
                self.label,
                weight=ft.FontWeight.BOLD,
                size=self.font_size_default,
            ),
            table_control,
        ]
        if action_controls:
            controls.append(
                ft.Row(
                    action_controls,
                    alignment=ft.MainAxisAlignment.CENTER,
                    spacing=4,
                    tight=True,
                )
            )
        super().__init__(
            controls=controls,
            spacing=spacing,
            horizontal_alignment=ft.CrossAxisAlignment.CENTER,
        )

    def _notify_changed(self) -> None:
        """Invoke the change handler and update the control if mounted."""
        if self.on_change_handler is not None:
            self.on_change_handler(None)
        if self.page:
            self.update()
