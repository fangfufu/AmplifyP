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

"""Header component for the primers list."""

from __future__ import annotations

from collections.abc import Callable

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings


class PrimerHeader(ft.Container):  # type: ignore[misc]
    """Header component for the primers list showing columns and resizing."""

    def __init__(
        self,
        settings: GUISettings,
        on_toggle_all: Callable[[ft.Event[ft.Checkbox]], None],
        on_divider_pan: Callable[[ft.DragUpdateEvent], None],
        on_divider_pan_end: Callable[[ft.DragEndEvent], None],
        name_column_width: float,
        on_add_primer: Callable[[ft.Event | None], None],
        on_delete_primer: Callable[[ft.Event | None], None],
        on_move_primer_up: Callable[[ft.Event | None], None],
        on_move_primer_down: Callable[[ft.Event | None], None],
    ) -> None:
        """Initialise the PrimerHeader."""
        self.settings = settings
        self.all_primers_checkbox = ft.Checkbox(
            value=None,
            tristate=True,
            on_change=on_toggle_all,
        )
        show_temp = self.settings.get("show_primer_temperature", False)
        self.add_button = ft.IconButton(
            icon=ft.Icons.ADD_CIRCLE_OUTLINE,
            icon_size=16,
            width=24,
            height=24,
            padding=0,
            tooltip="Add Primer Below",
            disabled=True,
            on_click=on_add_primer,
        )
        self.delete_button = ft.IconButton(
            icon=ft.Icons.DELETE_OUTLINE,
            icon_size=16,
            width=24,
            height=24,
            padding=0,
            tooltip="Delete Primer",
            disabled=True,
            on_click=on_delete_primer,
        )
        self.up_button = ft.IconButton(
            icon=ft.Icons.ARROW_UPWARD,
            icon_size=16,
            width=24,
            height=24,
            padding=0,
            tooltip="Move Up",
            disabled=True,
            on_click=on_move_primer_up,
        )
        self.down_button = ft.IconButton(
            icon=ft.Icons.ARROW_DOWNWARD,
            icon_size=16,
            width=24,
            height=24,
            padding=0,
            tooltip="Move Down",
            disabled=True,
            on_click=on_move_primer_down,
        )
        self.reorder_controls = ft.Row(
            [
                self.add_button,
                self.delete_button,
                self.up_button,
                self.down_button,
            ],
            spacing=2,
            alignment=ft.MainAxisAlignment.CENTER,
        )
        self.tm_divider = ft.Container(
            width=4,
            bgcolor=GUIColours.DIVIDER_GREY,
            margin=0,
            height=36,
            visible=show_temp,
        )

        if show_temp:
            self.tm_header = ft.Container(
                content=self.reorder_controls,
                width=108,
                padding=ft.Padding(5, 0, 5, 0),
                alignment=ft.Alignment(0, 0),
                visible=True,
            )
            self.sequence_header = ft.Container(
                content=ft.Text(
                    "Sequence",
                    weight=ft.FontWeight.BOLD,
                    size=self.settings.get("font_size_small", 12),
                ),
                expand=True,
                padding=ft.Padding(5, 0, 0, 0),
                alignment=ft.Alignment(-1, 0),
            )
        else:
            self.tm_header = ft.Container(
                content=ft.Text(
                    "Tm",
                    weight=ft.FontWeight.BOLD,
                    size=self.settings.get("font_size_small", 12),
                ),
                width=50,
                padding=ft.Padding(5, 0, 0, 0),
                alignment=ft.Alignment(-1, 0),
                visible=False,
            )
            self.sequence_header = ft.Container(
                content=ft.Row(
                    [
                        ft.Text(
                            "Sequence",
                            weight=ft.FontWeight.BOLD,
                            size=self.settings.get("font_size_small", 12),
                        ),
                        self.reorder_controls,
                    ],
                    alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                ),
                expand=True,
                padding=ft.Padding(5, 0, 5, 0),
                alignment=ft.Alignment(-1, 0),
            )

        controls = [
            ft.Container(width=25),
            ft.Container(
                content=self.all_primers_checkbox,
                width=30,
                alignment=ft.Alignment(0, 0),
            ),
            ft.Container(
                width=4,
                bgcolor=GUIColours.DIVIDER_GREY,
                margin=0,
                height=36,
            ),
            ft.Container(
                content=ft.Text(
                    "Name",
                    weight=ft.FontWeight.BOLD,
                    size=self.settings.get("font_size_small", 12),
                ),
                width=name_column_width,
                padding=ft.Padding(5, 0, 0, 0),
                alignment=ft.Alignment(-1, 0),
            ),
            ft.GestureDetector(
                on_pan_update=on_divider_pan,
                on_pan_end=on_divider_pan_end,
                content=ft.Container(
                    width=4,
                    bgcolor=GUIColours.DIVIDER_GREY,
                    margin=0,
                    height=36,
                ),
                mouse_cursor=ft.MouseCursor.RESIZE_LEFT_RIGHT,
            ),
            self.sequence_header,
        ]
        if show_temp:
            controls.append(self.tm_header)
        self.header_row = ft.Row(
            controls,
            alignment=ft.MainAxisAlignment.START,
            height=36,
            spacing=0,
        )
        super().__init__(
            content=self.header_row,
            padding=0,
            border=ft.Border(bottom=ft.BorderSide(4, GUIColours.DIVIDER_GREY)),
            height=40,
        )
