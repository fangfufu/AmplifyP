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

"""InputTile expansion tile component for settings view."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.gui_helpers import BorderedCheckbox


class InputTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Input settings (Primer List and Tm)."""

    def __init__(
        self,
        settings: GUISettings,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.Event | None], None],
        header_size: int,
        font_size_default: int,
    ) -> None:
        """Initialise the InputTile."""
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        # --- Primer List settings ---
        self.set_primer_info_panel_position = ft.Dropdown(
            label="Primer Info Panel Position",
            options=[
                ft.dropdown.Option(key="top", text="Top"),
                ft.dropdown.Option(key="bottom", text="Bottom"),
            ],
            width=700,
            on_select=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        self.fixed_height_primer_info_checkbox = BorderedCheckbox(
            label="Fixed height primer info box",
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
        self.settings_map["primer_info_panel_fixed_height"] = (
            self.fixed_height_primer_info_checkbox
        )
        self.settings_map["show_primer_temperature"] = (
            self.set_show_primer_temperature
        )
        self.settings_map["auto_activate_new_valid_primer"] = (
            self.auto_activate_new_valid_primer
        )
        self.settings_map["tm_colour_scheme"] = self.set_tm_colour_scheme

        # --- Primer Melting Temperature (Tm) settings ---
        self.set_tm_dna_conc = ft.TextField(
            label="DNA Conc (nM)",
            value="50.0",
            on_change=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )
        self.set_tm_dnap_conc = ft.TextField(
            label="DNA Pol Conc",
            value="0.0",
            on_change=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )
        self.set_tm_mono_salt = ft.TextField(
            label="Monovalent Salt Conc (mM)",
            value="50.0",
            on_change=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )
        self.set_tm_div_salt = ft.TextField(
            label="Divalent Salt Conc (mM)",
            value="1.5",
            on_change=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )
        self.set_tm_dNTP_conc = ft.TextField(
            label="dNTP Conc (mM)",
            value="0.0",
            on_change=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )
        self.set_tm_method = ft.Dropdown(
            label="Tm Calculation Method",
            options=[
                ft.dropdown.Option("SantaLucia 1998 / Owczarzy 2008 (Default)"),
                ft.dropdown.Option("Lander / Amplify 4"),
            ],
            expand=True,
            on_select=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )
        self.settings_map["tm_dna_conc"] = self.set_tm_dna_conc
        self.settings_map["tm_dnap_conc"] = self.set_tm_dnap_conc
        self.settings_map["tm_mono_salt"] = self.set_tm_mono_salt
        self.settings_map["tm_div_salt"] = self.set_tm_div_salt
        self.settings_map["tm_dNTP_conc"] = self.set_tm_dNTP_conc
        self.settings_map["tm_method"] = self.set_tm_method

        super().__init__(
            title=ft.Text(
                "Input",
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
                                        ft.Text(
                                            "Primer List",
                                            weight=ft.FontWeight.BOLD,
                                            size=font_size_default + 2,
                                        ),
                                        self.set_primer_info_panel_position,
                                        self.fixed_height_primer_info_checkbox,
                                        self.auto_activate_new_valid_primer,
                                        self.set_show_primer_temperature,
                                        ft.Row([self.set_tm_colour_scheme]),
                                        ft.Divider(),
                                        ft.Text(
                                            "Primer Melting Temperature (Tm)",
                                            weight=ft.FontWeight.BOLD,
                                            size=font_size_default + 2,
                                        ),
                                        ft.Row([self.set_tm_method]),
                                        self.set_tm_dna_conc,
                                        self.set_tm_dnap_conc,
                                        self.set_tm_mono_salt,
                                        self.set_tm_div_salt,
                                        self.set_tm_dNTP_conc,
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=700,
                            ),
                        ],
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )
