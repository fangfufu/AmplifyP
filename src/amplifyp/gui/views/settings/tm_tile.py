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

"""TmTile component for Flet settings view."""

from typing import Any

import flet as ft


class TmTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Primer Melting Temperature (Tm) settings."""

    def __init__(
        self,
        settings: Any,
        settings_map: dict[str, Any],
        on_change_handler: Any,
        header_size: int,
    ) -> None:
        """Initialise the TmTile.

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

        self.set_tm_dna_conc = ft.TextField(
            label="DNA Conc (nM)",
            value="50.0",
            on_change=self.on_change_handler,
        )
        self.set_tm_dnap_conc = ft.TextField(
            label="DNA Pol Conc",
            value="0.0",
            on_change=self.on_change_handler,
        )
        self.set_tm_mono_salt = ft.TextField(
            label="Monovalent Salt Conc (mM)",
            value="50.0",
            on_change=self.on_change_handler,
        )
        self.set_tm_div_salt = ft.TextField(
            label="Divalent Salt Conc (mM)",
            value="1.5",
            on_change=self.on_change_handler,
        )
        self.set_tm_dNTP_conc = ft.TextField(
            label="dNTP Conc (mM)",
            value="0.0",
            on_change=self.on_change_handler,
        )
        self.set_tm_method = ft.Dropdown(
            label="Tm Calculation Method",
            options=[
                ft.dropdown.Option("SantaLucia 1998 / Owczarzy 2008 (Default)"),
                ft.dropdown.Option("Lander / Amplify 4"),
            ],
            expand=True,
            on_select=self.on_change_handler,
        )
        from amplifyp.gui.util import BorderedCheckbox

        self.set_show_primer_temperature = BorderedCheckbox(
            label="Show Primer Temperature in List",
            on_change=self.on_change_handler,
        )

        self.settings_map["tm_dna_conc"] = self.set_tm_dna_conc
        self.settings_map["tm_dnap_conc"] = self.set_tm_dnap_conc
        self.settings_map["tm_mono_salt"] = self.set_tm_mono_salt
        self.settings_map["tm_div_salt"] = self.set_tm_div_salt
        self.settings_map["tm_dNTP_conc"] = self.set_tm_dNTP_conc
        self.settings_map["tm_method"] = self.set_tm_method
        self.settings_map["show_primer_temperature"] = (
            self.set_show_primer_temperature
        )

        super().__init__(
            title=ft.Text(
                "Primer Melting Temperature (Tm) Settings",
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
                                    ft.Container(expand=True),
                                    ft.Container(
                                        content=ft.Column(
                                            [
                                                ft.Row([self.set_tm_method]),
                                                self.set_tm_dna_conc,
                                                self.set_tm_dnap_conc,
                                                self.set_tm_mono_salt,
                                                self.set_tm_div_salt,
                                                self.set_tm_dNTP_conc,
                                                self.set_show_primer_temperature,
                                            ],
                                            spacing=15,
                                            horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                        ),
                                        expand=2,
                                    ),
                                    ft.Container(expand=True),
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
