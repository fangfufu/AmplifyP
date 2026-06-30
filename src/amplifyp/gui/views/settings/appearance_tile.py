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

"""AppearanceTile component for Flet settings view."""

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings


class AppearanceTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Appearance settings."""

    def __init__(
        self,
        settings: GUISettings,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.Event | None], None],
        header_size: int,
    ) -> None:
        """Initialise the AppearanceTile.

        Args:
            settings: The settings object.
            settings_map: A dictionary mapping setting keys to UI
                components for population and retrieval.
            on_change_handler: The handler to call when a setting changes.
            header_size: The size of the header text.
        """
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
            width=350,
            on_select=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        self.set_colour_scheme = ft.Dropdown(
            label="Colour Scheme",
            options=[
                ft.dropdown.Option("Light"),
                ft.dropdown.Option("Dark"),
                ft.dropdown.Option("Light (Colour Deficient Friendly)"),
                ft.dropdown.Option("Dark (Colour Deficient Friendly)"),
                ft.dropdown.Option("System"),
                ft.dropdown.Option("System (Colour Deficient Friendly)"),
            ],
            width=350,
            on_select=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        self._dummy_colour_deficient = ft.Checkbox(visible=False)

        self.settings_map["font_family"] = self.set_font_family
        self.settings_map["colour_deficient"] = self._dummy_colour_deficient

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
                                        self.set_colour_scheme,
                                        self._dummy_colour_deficient,
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=350,
                            ),
                        ],
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )

    @property
    def set_colour_deficient(self) -> ft.Checkbox:
        """Get colour deficient checkbox for backwards compatibility."""
        return self._dummy_colour_deficient

    def _on_colour_scheme_change(self, e: ft.ControlEvent) -> None:
        """Handle colour scheme dropdown change.

        Syncs the selected colour scheme to settings and triggers the
        change handler.

        Args:
            e: The Flet control event triggered by the dropdown change.
        """
        self.sync_colour_scheme_to_settings()
        self.on_change_handler(e)

    def update_colour_scheme_dropdown(self) -> None:
        """Update the colour scheme dropdown value based on settings.

        Reads the current dark_mode and colour_deficient settings and
        sets the dropdown to the matching display option.
        """
        dark = self.settings.get("dark_mode", False)
        deficient = bool(self.settings.get("colour_deficient", False))
        if str(dark).lower() == "system":
            if deficient:
                self.set_colour_scheme.value = (
                    "System (Colour Deficient Friendly)"
                )
            else:
                self.set_colour_scheme.value = "System"
        elif bool(dark):
            if deficient:
                self.set_colour_scheme.value = (
                    "Dark (Colour Deficient Friendly)"
                )
            else:
                self.set_colour_scheme.value = "Dark"
        else:
            if deficient:
                self.set_colour_scheme.value = (
                    "Light (Colour Deficient Friendly)"
                )
            else:
                self.set_colour_scheme.value = "Light"

    def sync_colour_scheme_to_settings(self) -> None:
        """Sync the colour scheme dropdown selection back to settings.

        Parses the dropdown value and updates dark_mode and
        colour_deficient settings accordingly, including the hidden
        checkbox state.
        """
        val = self.set_colour_scheme.value
        if val == "Dark":
            self.settings["dark_mode"] = True
            self.settings["colour_deficient"] = False
            self._dummy_colour_deficient.value = False
        elif val == "Dark (Colour Deficient Friendly)":
            self.settings["dark_mode"] = True
            self.settings["colour_deficient"] = True
            self._dummy_colour_deficient.value = True
        elif val == "Light (Colour Deficient Friendly)":
            self.settings["dark_mode"] = False
            self.settings["colour_deficient"] = True
            self._dummy_colour_deficient.value = True
        elif val == "System":
            self.settings["dark_mode"] = "system"
            self.settings["colour_deficient"] = False
            self._dummy_colour_deficient.value = False
        elif val == "System (Colour Deficient Friendly)":
            self.settings["dark_mode"] = "system"
            self.settings["colour_deficient"] = True
            self._dummy_colour_deficient.value = True
        else:  # Light
            self.settings["dark_mode"] = False
            self.settings["colour_deficient"] = False
            self._dummy_colour_deficient.value = False
