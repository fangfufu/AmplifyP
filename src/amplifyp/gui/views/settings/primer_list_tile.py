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

"""PrimerListTile component for Flet settings view."""

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.ui import BorderedCheckbox


class PrimerListTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Primer List settings."""

    def __init__(
        self,
        settings: GUISettings,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.Event | None], None],
        header_size: int,
    ) -> None:
        """Initialise the PrimerListTile.

        Args:
            settings: The settings object.
            settings_map: A dictionary mapping setting keys to UI
                components for population and retrieval.
            on_change_handler: The handler to call when a setting changes.
            header_size: The size of the expansion tile header text.
        """
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        self.set_primer_info_panel_position = ft.Dropdown(
            label="Primer Info Panel Position",
            options=[
                ft.dropdown.Option(key="top", text="Top"),
                ft.dropdown.Option(key="bottom", text="Bottom"),
            ],
            width=500,
            on_select=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        self.ignore_inactive_name_dup_checkbox = BorderedCheckbox(
            label=("Ignore inactive primers when checking for duplicate names"),
            on_change=self.on_change_handler,
        )

        self.ignore_inactive_seq_dup_checkbox = BorderedCheckbox(
            label=(
                "Ignore inactive primers when checking for duplicate sequences"
            ),
            on_change=self.on_change_handler,
        )

        self.set_show_primer_temperature = BorderedCheckbox(
            label="Show primer temperature column",
            on_change=self.on_change_handler,
        )

        self.auto_activate_new_valid_primer = BorderedCheckbox(
            label="Auto-activate new valid primer",
            on_change=self.on_change_handler,
        )

        self.set_tm_colour_scheme = ft.Dropdown(
            label="Primer temperature colour scheme",
            options=[
                ft.dropdown.Option("None"),
                ft.dropdown.Option("Cool-Warm"),
                ft.dropdown.Option("Traffic Light"),
            ],
            expand=True,
            on_select=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        self.settings_map["primer_info_panel_position"] = (
            self.set_primer_info_panel_position
        )
        self.settings_map["ignore_inactive_name_dup_warn"] = (
            self.ignore_inactive_name_dup_checkbox
        )
        self.settings_map["ignore_inactive_seq_dup_warn"] = (
            self.ignore_inactive_seq_dup_checkbox
        )
        self.settings_map["show_primer_temperature"] = (
            self.set_show_primer_temperature
        )
        self.settings_map["auto_activate_new_valid_primer"] = (
            self.auto_activate_new_valid_primer
        )
        self.settings_map["tm_colour_scheme"] = self.set_tm_colour_scheme

        super().__init__(
            title=ft.Text(
                "Primer List Settings",
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
                                        self.set_primer_info_panel_position,
                                        self.ignore_inactive_name_dup_checkbox,
                                        self.ignore_inactive_seq_dup_checkbox,
                                        self.auto_activate_new_valid_primer,
                                        self.set_show_primer_temperature,
                                        ft.Row([self.set_tm_colour_scheme]),
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=500,
                            ),
                        ],
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )
