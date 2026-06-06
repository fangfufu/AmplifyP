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

"""Settings View for the Flet application."""

from typing import Any

import flet as ft

from amplifyp.gui.state import GUIColors, GUIState
from amplifyp.settings import ReplicationSettings


class SettingsView(ft.ListView):  # type: ignore[misc]
    """Settings view for configuring PCR replication settings."""

    def __init__(
        self,
        page: ft.Page,
        state: GUIState | None = None,
        on_change: Any | None = None,
        on_apply: Any | None = None,
        on_reset: Any | None = None,
    ) -> None:
        """Initialize the SettingsView."""
        super().__init__(expand=True, spacing=20, padding=10)
        self.app_page = page
        self.state = state if state is not None else GUIState()
        self.on_change = on_change
        self.on_apply = on_apply
        self.on_reset = on_reset

        # Settings State
        # Replication Settings
        self.set_primability_cutoff = ft.TextField(
            label="Primability Cutoff",
            value="0.8",
            on_change=self.on_change_handler,
        )
        self.set_stability_cutoff = ft.TextField(
            label="Stability Cutoff",
            value="0.4",
            on_change=self.on_change_handler,
        )
        self.set_amp4_compat = ft.Checkbox(
            label="Amplify4 Compatibility Mode",
            value=False,
            on_change=self.on_change_handler,
        )

        # Primer Melting Temperature (Tm) Settings
        self.set_tm_dna_conc = ft.TextField(
            label="DNA Conc (nM)",
            value="50.0",
            on_change=self.on_change_handler,
        )
        self.set_tm_dnap_conc = ft.TextField(
            label="DNA Pol Conc", value="0.0", on_change=self.on_change_handler
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
                ft.dropdown.Option("Amplify P (Default)"),
                ft.dropdown.Option("Amplify 4"),
            ],
        )
        self.set_tm_method.on_change = self.on_change_handler

        # Amplify4 TM Settings
        self.set_amp4tm_dna_conc = ft.TextField(
            label="DNA Conc (nM)",
            value="50.0",
            on_change=self.on_change_handler,
        )
        self.set_amp4tm_mono_salt = ft.TextField(
            label="Monovalent Salt Conc (mM)",
            value="50.0",
            on_change=self.on_change_handler,
        )

        # Primer Dimer Settings
        self.set_pd_min_overlap = ft.TextField(
            label="Min Overlap", value="3", on_change=self.on_change_handler
        )
        self.set_pd_threshold = ft.TextField(
            label="Threshold", value="60.0", on_change=self.on_change_handler
        )

        # Appearance Settings
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

        self.settings_map = {
            "primability_cutoff": self.set_primability_cutoff,
            "stability_cutoff": self.set_stability_cutoff,
            "amp4_compat": self.set_amp4_compat,
            "tm_dna_conc": self.set_tm_dna_conc,
            "tm_dnap_conc": self.set_tm_dnap_conc,
            "tm_mono_salt": self.set_tm_mono_salt,
            "tm_div_salt": self.set_tm_div_salt,
            "tm_dNTP_conc": self.set_tm_dNTP_conc,
            "tm_method": self.set_tm_method,
            "amp4tm_dna_conc": self.set_amp4tm_dna_conc,
            "amp4tm_mono_salt": self.set_amp4tm_mono_salt,
            "pd_min_overlap": self.set_pd_min_overlap,
            "pd_threshold": self.set_pd_threshold,
            "font_family": self.set_font_family,
            "color_deficient": self.set_color_deficient,
        }

        # Dynamic Base Pair Scores settings mapping
        from amplifyp.dna import Nucleotides

        for r_char in Nucleotides.PRIMER:
            for c_char in Nucleotides.TEMPLATE:
                if c_char == Nucleotides.GAP:
                    continue
                key = f"bp_score_{r_char}_{c_char}"
                self.settings_map[key] = ft.TextField(
                    value="0",
                    on_change=self.on_change_handler,
                    text_align=ft.TextAlign.CENTER,
                    dense=True,
                    width=60,
                    height=36,
                    content_padding=4,
                    text_style=ft.TextStyle(color=ft.Colors.BLACK, size=14),
                )

        # Build Base Pair Scores Styled Table (matching primer table)
        header_controls = [
            ft.Container(
                content=ft.Text(
                    "Primer / Template",
                    weight=ft.FontWeight.BOLD,
                    size=11,
                ),
                width=120,
                alignment=ft.Alignment(-1, 0),
                padding=ft.Padding(10, 0, 0, 0),
            )
        ]

        for c_char in Nucleotides.TEMPLATE:
            if c_char == Nucleotides.GAP:
                continue
            header_controls.append(
                ft.Container(
                    width=1,
                    bgcolor=GUIColors.DIVIDER_GREY,
                    height=36,
                )
            )
            header_controls.append(
                ft.Container(
                    content=ft.Text(c_char, weight=ft.FontWeight.BOLD, size=15),
                    width=60,
                    alignment=ft.Alignment(0, 0),
                )
            )

        header_row = ft.Container(
            content=ft.Row(
                header_controls,
                spacing=0,
                alignment=ft.MainAxisAlignment.START,
            ),
            border=ft.Border(bottom=ft.BorderSide(2, GUIColors.DIVIDER_GREY)),
            height=36,
        )

        row_controls = [header_row]
        for r_char in Nucleotides.PRIMER:
            cells = [
                ft.Container(
                    content=ft.Text(r_char, weight=ft.FontWeight.BOLD, size=15),
                    width=120,
                    alignment=ft.Alignment(-1, 0),
                    padding=ft.Padding(10, 0, 0, 0),
                )
            ]
            for c_char in Nucleotides.TEMPLATE:
                if c_char == Nucleotides.GAP:
                    continue
                key = f"bp_score_{r_char}_{c_char}"

                # Style TextField to match the borderless style in primer table
                self.settings_map[key].border = ft.InputBorder.NONE
                self.settings_map[key].width = 60
                self.settings_map[key].height = 30
                self.settings_map[key].content_padding = 0

                cells.append(
                    ft.Container(
                        width=1,
                        bgcolor=GUIColors.DIVIDER_GREY,
                        height=30,
                    )
                )
                cells.append(
                    ft.Container(
                        content=self.settings_map[key],
                        width=60,
                        alignment=ft.Alignment(0, 0),
                    )
                )

            row_controls.append(
                ft.Container(
                    content=ft.Row(
                        cells,
                        spacing=0,
                        alignment=ft.MainAxisAlignment.START,
                    ),
                    border=ft.Border(
                        bottom=ft.BorderSide(1, GUIColors.DIVIDER_GREY)
                    ),
                    height=30,
                )
            )

        bp_table_container = ft.Column(
            [
                ft.Text(
                    "Base Pair Weights", weight=ft.FontWeight.BOLD, size=14
                ),
                ft.Row(
                    [
                        ft.Container(
                            content=ft.Column(
                                row_controls,
                                spacing=0,
                            ),
                            border=ft.Border.all(1, GUIColors.OUTLINE),
                            border_radius=5,
                            padding=0,
                        )
                    ],
                    scroll=ft.ScrollMode.ADAPTIVE,
                ),
            ],
            spacing=10,
        )

        header_size = self.state.settings.get("font_size_header", 18)
        self.controls = [
            ft.ExpansionTile(
                title=ft.Text(
                    "Origin of Replication Settings",
                    weight=ft.FontWeight.BOLD,
                    size=header_size,
                ),
                expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
                controls=[
                    ft.Container(
                        content=ft.Column(
                            [
                                self.set_primability_cutoff,
                                self.set_stability_cutoff,
                                self.set_amp4_compat,
                                ft.Divider(),
                                bp_table_container,
                            ],
                            spacing=15,
                            horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                        ),
                        padding=ft.Padding(0, 20, 0, 10),
                    )
                ],
            ),
            ft.ExpansionTile(
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
                                self.set_tm_method,
                                ft.Row(
                                    [
                                        ft.Container(
                                            content=ft.Column(
                                                [
                                                    ft.Text(
                                                        "Amplify P",
                                                        weight=ft.FontWeight.BOLD,
                                                    ),
                                                    self.set_tm_dna_conc,
                                                    self.set_tm_dnap_conc,
                                                    self.set_tm_mono_salt,
                                                    self.set_tm_div_salt,
                                                    self.set_tm_dNTP_conc,
                                                ],
                                                spacing=15,
                                                horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                            ),
                                            expand=True,
                                        ),
                                        ft.VerticalDivider(),
                                        ft.Container(
                                            content=ft.Column(
                                                [
                                                    ft.Text(
                                                        "Amplify 4",
                                                        weight=ft.FontWeight.BOLD,
                                                    ),
                                                    self.set_amp4tm_dna_conc,
                                                    self.set_amp4tm_mono_salt,
                                                ],
                                                spacing=15,
                                                horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                            ),
                                            expand=True,
                                        ),
                                    ],
                                    vertical_alignment=ft.CrossAxisAlignment.START,
                                    expand=True,
                                ),
                            ],
                            spacing=15,
                            horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                        ),
                        padding=ft.Padding(0, 20, 0, 10),
                    )
                ],
            ),
            ft.ExpansionTile(
                title=ft.Text(
                    "Primer Dimer Settings",
                    weight=ft.FontWeight.BOLD,
                    size=header_size,
                ),
                expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
                controls=[
                    ft.Container(
                        content=ft.Column(
                            [
                                self.set_pd_min_overlap,
                                self.set_pd_threshold,
                            ],
                            spacing=15,
                            horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                        ),
                        padding=ft.Padding(0, 20, 0, 10),
                    )
                ],
            ),
            ft.ExpansionTile(
                title=ft.Text(
                    "Appearance Settings",
                    weight=ft.FontWeight.BOLD,
                    size=header_size,
                ),
                expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
                controls=[
                    ft.Container(
                        content=ft.Column(
                            [
                                self.set_font_family,
                                self.set_color_deficient,
                            ],
                            spacing=15,
                            horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                        ),
                        padding=ft.Padding(0, 20, 0, 10),
                    )
                ],
            ),
            ft.Divider(),
            ft.Row(
                [
                    ft.FilledButton(
                        "Apply",
                        icon=ft.Icons.DONE,
                        on_click=self.on_apply_handler,
                    ),
                    ft.OutlinedButton(
                        "Reset to Default",
                        icon=ft.Icons.RESTORE,
                        on_click=self.on_reset_handler,
                    ),
                ],
                spacing=10,
            ),
        ]

        # Sync initial UI state
        self.update_ui()

    def sync_to_state(self) -> None:
        """Sync current settings UI controls back to the central state."""
        for k, v in self.settings_map.items():
            self.state.settings[k] = v.value

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state."""
        for k, v in self.state.settings.items():
            if k in self.settings_map:
                if isinstance(self.settings_map[k], ft.Checkbox):
                    self.settings_map[k].value = bool(v)
                else:
                    self.settings_map[k].value = str(v)

    def on_change_handler(self, e: ft.ControlEvent) -> None:
        """Handle change in settings fields."""
        self.sync_to_state()
        if self.state.has_run_pcr:
            self.state.results_outdated = True
        if self.on_change:
            self.on_change(e)

    def on_apply_handler(self, e: ft.ControlEvent) -> None:
        """Handle apply button click."""
        self.sync_to_state()
        if self.on_apply:
            self.on_apply(e)

    def on_reset_handler(self, e: ft.ControlEvent) -> None:
        """Handle reset to default button click."""
        from amplifyp.dna import Nucleotides
        from amplifyp.settings import (
            DEFAULT_BASE_PAIR_WEIGHTS,
            DEFAULT_PRIMABILITY_CUTOFF,
            DEFAULT_PRIMER_DIMER_OVERLAP,
            DEFAULT_PRIMER_DIMER_THRESHOLD,
            DEFAULT_STABILITY_CUTOFF,
            GLOBAL_AMPLIFY4_TM_SETTINGS,
            GLOBAL_TM_SETTINGS,
        )

        reset_dict = {
            "primability_cutoff": str(DEFAULT_PRIMABILITY_CUTOFF),
            "stability_cutoff": str(DEFAULT_STABILITY_CUTOFF),
            "amp4_compat": False,
            "tm_method": "Amplify P (Default)",
            "tm_dna_conc": str(GLOBAL_TM_SETTINGS.dna_conc),
            "tm_dnap_conc": str(GLOBAL_TM_SETTINGS.dnap_conc),
            "tm_mono_salt": str(GLOBAL_TM_SETTINGS.monovalent_salt_conc),
            "tm_div_salt": str(GLOBAL_TM_SETTINGS.divalent_salt_conc),
            "tm_dNTP_conc": str(GLOBAL_TM_SETTINGS.dnTP_conc),
            "amp4tm_dna_conc": str(GLOBAL_AMPLIFY4_TM_SETTINGS.dna_conc),
            "amp4tm_mono_salt": str(
                GLOBAL_AMPLIFY4_TM_SETTINGS.monovalent_salt_conc
            ),
            "pd_min_overlap": str(DEFAULT_PRIMER_DIMER_OVERLAP),
            "pd_threshold": str(DEFAULT_PRIMER_DIMER_THRESHOLD),
            "font_family": "Roboto Mono",
            "color_deficient": False,
        }

        for r_char in Nucleotides.PRIMER:
            for c_char in Nucleotides.TEMPLATE:
                if c_char == Nucleotides.GAP:
                    continue
                key = f"bp_score_{r_char}_{c_char}"
                reset_dict[key] = str(DEFAULT_BASE_PAIR_WEIGHTS[r_char, c_char])

        self.state.settings = reset_dict
        self.update_ui()
        self.app_page.update()
        if self.on_reset:
            self.on_reset(e)

    def get_replication_settings(self) -> ReplicationSettings:
        """Get the current settings as a ReplicationSettings object."""
        self.sync_to_state()
        return self.state.get_replication_settings()

    def get_state(self) -> dict[str, Any]:
        """Get the current state for serialization."""
        self.sync_to_state()
        return dict(self.state.settings)

    def set_state(self, state: dict[str, Any]) -> None:
        """Set the current state from deserialized data."""
        self.state.from_dict(state)
        self.update_ui()
        self.app_page.update()
