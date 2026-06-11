# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""AppearanceTile component for Flet settings view."""

from typing import Any

import flet as ft


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

        self.set_color_scheme = ft.Dropdown(
            label="Colour Scheme",
            options=[
                ft.dropdown.Option("Light"),
                ft.dropdown.Option("Dark"),
                ft.dropdown.Option("Light (Colour Deficient Friendly)"),
                ft.dropdown.Option("Dark (Colour Deficient Friendly)"),
            ],
        )
        self.set_color_scheme.on_change = self._on_color_scheme_change

        self._dummy_color_deficient = ft.Checkbox(visible=False)

        self.settings_map["font_family"] = self.set_font_family
        self.settings_map["color_deficient"] = self._dummy_color_deficient

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
                                        self.set_color_scheme,
                                        self._dummy_color_deficient,
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=450,
                            ),
                        ],
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )

    @property
    def set_color_deficient(self) -> ft.Checkbox:
        """Get color deficient checkbox for backwards compatibility."""
        return self._dummy_color_deficient

    def _on_color_scheme_change(self, e: ft.ControlEvent) -> None:
        """Handle color scheme dropdown change."""
        self.sync_color_scheme_to_settings()
        if self.on_change_handler:
            self.on_change_handler(e)

    def update_color_scheme_dropdown(self) -> None:
        """Update the color scheme dropdown value based on settings."""
        dark = bool(self.settings.get("dark_mode", False))
        deficient = bool(self.settings.get("color_deficient", False))
        if dark:
            if deficient:
                self.set_color_scheme.value = "Dark (Colour Deficient Friendly)"
            else:
                self.set_color_scheme.value = "Dark"
        else:
            if deficient:
                self.set_color_scheme.value = (
                    "Light (Colour Deficient Friendly)"
                )
            else:
                self.set_color_scheme.value = "Light"

    def sync_color_scheme_to_settings(self) -> None:
        """Sync the color scheme dropdown selection back to settings."""
        val = self.set_color_scheme.value
        if val == "Dark":
            self.settings["dark_mode"] = True
            self.settings["color_deficient"] = False
            self._dummy_color_deficient.value = False
        elif val == "Dark (Colour Deficient Friendly)":
            self.settings["dark_mode"] = True
            self.settings["color_deficient"] = True
            self._dummy_color_deficient.value = True
        elif val == "Light (Colour Deficient Friendly)":
            self.settings["dark_mode"] = False
            self.settings["color_deficient"] = True
            self._dummy_color_deficient.value = True
        else:  # Light
            self.settings["dark_mode"] = False
            self.settings["color_deficient"] = False
            self._dummy_color_deficient.value = False
