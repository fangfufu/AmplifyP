# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""TmTile component for Flet settings view."""

import flet as ft
from typing import Any


class TmTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Primer Melting Temperature (Tm) settings."""

    def __init__(
        self,
        settings: Any,
        settings_map: dict[str, Any],
        on_change_handler: Any,
        header_size: int,
    ) -> None:
        """Initialize the TmTile."""
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
        )
        self.set_tm_method.on_change = self.on_change_handler

        self.settings_map["tm_dna_conc"] = self.set_tm_dna_conc
        self.settings_map["tm_dnap_conc"] = self.set_tm_dnap_conc
        self.settings_map["tm_mono_salt"] = self.set_tm_mono_salt
        self.settings_map["tm_div_salt"] = self.set_tm_div_salt
        self.settings_map["tm_dNTP_conc"] = self.set_tm_dNTP_conc
        self.settings_map["tm_method"] = self.set_tm_method

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
