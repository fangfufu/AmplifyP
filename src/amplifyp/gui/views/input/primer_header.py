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

from amplifyp.gui.settings import GUIColors, GUISettings


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
        """Initialize the PrimerHeader."""
        self.settings = settings
        self.all_primers_checkbox = ft.Checkbox(
            value=None,
            tristate=True,
            on_change=on_toggle_all,
        )
        self.header_row = ft.Row(
            [
                ft.Container(
                    content=self.all_primers_checkbox,
                    width=55,
                    alignment=ft.Alignment(0, 0),
                ),
                ft.Container(
                    width=4,
                    bgcolor=GUIColors.DIVIDER_GREY,
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
                        bgcolor=GUIColors.DIVIDER_GREY,
                        margin=0,
                        height=36,
                    ),
                    mouse_cursor=ft.MouseCursor.RESIZE_LEFT_RIGHT,
                ),
                ft.Container(
                    content=ft.Text(
                        "Sequence",
                        weight=ft.FontWeight.BOLD,
                        size=self.settings.get("font_size_small", 12),
                    ),
                    expand=True,
                    padding=ft.Padding(5, 0, 0, 0),
                    alignment=ft.Alignment(-1, 0),
                ),
            ],
            alignment=ft.MainAxisAlignment.START,
            height=36,
            spacing=0,
        )
        super().__init__(
            content=self.header_row,
            padding=0,
            border=ft.Border(bottom=ft.BorderSide(4, GUIColors.DIVIDER_GREY)),
            height=40,
        )
