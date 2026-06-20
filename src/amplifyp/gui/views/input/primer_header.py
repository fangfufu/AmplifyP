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

"""Header component for the primers list showing columns and resizing."""

from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings


class PrimerHeader(ft.Container):  # type: ignore[misc]
    """Header component for the primers list showing columns and resizing."""

    def __init__(
        self,
        settings: GUISettings,
        on_toggle_all: Any,
        on_divider_pan: Any,
        on_divider_pan_end: Any,
        name_column_width: float,
    ) -> None:
        """Initialise the PrimerHeader.

        Args:
            settings: Application GUI settings instance.
            on_toggle_all: Callback for the tri-state "all primers" checkbox.
            on_divider_pan: Callback for dragging the name/sequence divider.
            on_divider_pan_end: Callback for ending the divider drag.
            name_column_width: Width of the name column in pixels.
        """
        self.settings = settings
        self.all_primers_checkbox = ft.Checkbox(
            value=None,
            tristate=True,
            on_change=on_toggle_all,
        )
        show_temp = self.settings.get("show_primer_temperature", False)
        self.tm_header = ft.Container(
            content=ft.Text(
                "Tm",
                weight=ft.FontWeight.BOLD,
                size=self.settings.get("font_size_small", 12),
            ),
            width=50,
            padding=ft.Padding(5, 0, 0, 0),
            alignment=ft.Alignment(-1, 0),
            visible=show_temp,
        )
        self.tm_divider = ft.Container(
            width=4,
            bgcolor=GUIColours.DIVIDER_GREY,
            margin=0,
            height=36,
            visible=show_temp,
        )
        controls = [
            ft.Container(
                content=self.all_primers_checkbox,
                width=55,
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
        ]
        if show_temp:
            controls.extend([self.tm_header, self.tm_divider])
        controls.append(
            ft.Container(
                content=ft.Text(
                    "Sequence",
                    weight=ft.FontWeight.BOLD,
                    size=self.settings.get("font_size_small", 12),
                ),
                expand=True,
                padding=ft.Padding(5, 0, 0, 0),
                alignment=ft.Alignment(-1, 0),
            )
        )
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
