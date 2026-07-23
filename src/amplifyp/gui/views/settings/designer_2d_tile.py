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

"""Designer2DTile expansion tile component for settings view."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings


class Designer2DTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Designer 2D settings."""

    def __init__(
        self,
        settings: GUISettings,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.ControlEvent], None],
        header_size: int,
    ) -> None:
        """Initialise the Designer2DTile."""
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        self.set_designer_2d_colour_scheme = ft.Dropdown(
            label="2D Results Grid colour map scheme",
            options=[
                ft.dropdown.Option("None"),
                ft.dropdown.Option("Cool-Warm"),
                ft.dropdown.Option("Traffic Light"),
                ft.dropdown.Option("Blue-Orange"),
            ],
            expand=True,
            on_select=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        self.settings_map["designer_2d_colour_scheme"] = (
            self.set_designer_2d_colour_scheme
        )

        super().__init__(
            title=ft.Text(
                "Designer 2D Settings",
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Column(
                        [
                            ft.Row(
                                [
                                    ft.Container(
                                        content=ft.Column(
                                            [
                                                ft.Row(
                                                    [
                                                        self.set_designer_2d_colour_scheme
                                                    ]
                                                ),
                                            ],
                                            spacing=15,
                                            horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                        ),
                                        width=500,
                                    ),
                                ],
                                alignment=ft.MainAxisAlignment.CENTER,
                            ),
                        ],
                        spacing=15,
                        horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )
