# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""AppearanceTile component for Flet settings view."""

import flet as ft
from typing import Any


class AppearanceTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Appearance settings."""

    def __init__(
        self,
        settings: Any,
        settings_map: dict[str, Any],
        on_change_handler: Any,
        header_size: int,
    ) -> None:
        """Initialize the AppearanceTile."""
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        self.set_font_family = ft.Dropdown(
            label="Font Family",
            options=[
                ft.dropdown.Option("Roboto Mono"),
                ft.dropdown.Option("Courier New"),
                ft.dropdown.Option("Consolas"),
                ft.dropdown.Option("monospace"),
            ],
        )
        self.set_font_family.on_change = self.on_change_handler

        self.set_color_deficient = ft.Checkbox(
            label="Color Deficient Friendly Colour Scheme",
            value=False,
            on_change=self.on_change_handler,
        )

        self.settings_map["font_family"] = self.set_font_family
        self.settings_map["color_deficient"] = self.set_color_deficient

        super().__init__(
            title=ft.Text(
                "Appearance Settings",
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Row(
                        [
                            ft.Container(
                                content=ft.Column(
                                    [
                                        self.set_font_family,
                                        self.set_color_deficient,
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=300,
                            ),
                        ],
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )
