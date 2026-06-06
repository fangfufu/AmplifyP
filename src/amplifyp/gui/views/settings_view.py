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

from amplifyp.gui.settings import GUIColors, GUISettings
from amplifyp.settings import ReplicationSettings


class SettingsView(ft.ListView):  # type: ignore[misc]
    """Settings view for configuring PCR replication settings."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings | None = None,
        on_change: Any | None = None,
        on_apply: Any | None = None,
        on_reset: Any | None = None,
    ) -> None:
        """Initialize the SettingsView."""
        super().__init__(
            expand=True, spacing=20, padding=10, scroll=ft.ScrollMode.ALWAYS
        )
        self.app_page = page
        self.settings = settings if settings is not None else GUISettings()
        self.on_change = on_change
        self.on_apply = on_apply
        self.on_reset = on_reset

        self._init_basic_controls()

        font_size_default = self.settings.get("font_size_default", 14)
        font_size_micro = self.settings.get("font_size_micro", 10)
        font_size_table_header = self.settings.get("font_size_table_header", 15)

        self._init_dynamic_scores(font_size_default)

        bp_table_container = self._build_bp_table(
            font_size_default, font_size_micro, font_size_table_header
        )
        pd_table_container = self._build_pd_table(
            font_size_default, font_size_micro, font_size_table_header
        )

        header_size = self.settings.get("font_size_header", 18)
        self._build_layout(bp_table_container, pd_table_container, header_size)

        # Sync initial UI state
        self.update_ui()

    def _init_basic_controls(self) -> None:
        """Initialize basic control components and settings map."""
        # Replication Settings
        self.set_primability_cutoff = ft.TextField(
            label="Primability Cutoff",
            value="0.8",
            on_change=self._on_change_handler,
        )
        self.set_stability_cutoff = ft.TextField(
            label="Stability Cutoff",
            value="0.4",
            on_change=self._on_change_handler,
        )
        self.set_amp4_compat = ft.Checkbox(
            label="Amplify4 Compatibility Mode",
            value=False,
            on_change=self._on_change_handler,
        )

        # Primer Melting Temperature (Tm) Settings
        self.set_tm_dna_conc = ft.TextField(
            label="DNA Conc (nM)",
            value="50.0",
            on_change=self._on_change_handler,
        )
        self.set_tm_dnap_conc = ft.TextField(
            label="DNA Pol Conc", value="0.0", on_change=self._on_change_handler
        )
        self.set_tm_mono_salt = ft.TextField(
            label="Monovalent Salt Conc (mM)",
            value="50.0",
            on_change=self._on_change_handler,
        )
        self.set_tm_div_salt = ft.TextField(
            label="Divalent Salt Conc (mM)",
            value="1.5",
            on_change=self._on_change_handler,
        )
        self.set_tm_dNTP_conc = ft.TextField(
            label="dNTP Conc (mM)",
            value="0.0",
            on_change=self._on_change_handler,
        )
        self.set_tm_method = ft.Dropdown(
            label="Tm Calculation Method",
            options=[
                ft.dropdown.Option("SantaLucia 1998 / Owczarzy 2008 (Default)"),
                ft.dropdown.Option("Lander / Amplify 4"),
            ],
            expand=True,
        )
        self.set_tm_method.on_change = self._on_change_handler

        # Primer Dimer Settings
        self.set_pd_min_overlap = ft.TextField(
            label="Min Overlap", value="3", on_change=self._on_change_handler
        )
        self.set_pd_threshold = ft.TextField(
            label="Threshold", value="60.0", on_change=self._on_change_handler
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
        self.set_font_family.on_change = self._on_change_handler

        self.set_color_deficient = ft.Checkbox(
            label="Color Deficient Friendly Colour Scheme",
            value=False,
            on_change=self._on_change_handler,
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
            "pd_min_overlap": self.set_pd_min_overlap,
            "pd_threshold": self.set_pd_threshold,
            "font_family": self.set_font_family,
            "color_deficient": self.set_color_deficient,
        }

    def _init_dynamic_scores(self, font_size_default: int) -> None:
        """Initialize dynamic BP and PD score text fields in settings map."""
        from amplifyp.dna import Nucleotides

        # Dynamic Base Pair Scores settings mapping
        for r_char in Nucleotides.PRIMER:
            for c_char in Nucleotides.TEMPLATE:
                if c_char == Nucleotides.GAP:
                    continue
                key = f"bp_score_{r_char}_{c_char}"
                self.settings_map[key] = ft.TextField(
                    value="0",
                    on_change=self._on_change_handler,
                    text_align=ft.TextAlign.CENTER,
                    dense=True,
                    width=48,
                    height=36,
                    content_padding=4,
                    text_style=ft.TextStyle(
                        color=GUIColors.DIAGRAM_BLACK, size=font_size_default
                    ),
                )

        # Dynamic Primer Dimer Scores settings mapping
        for r_char in Nucleotides.PRIMER:
            for c_char in Nucleotides.PRIMER:
                key = f"pd_score_{r_char}_{c_char}"
                self.settings_map[key] = ft.TextField(
                    value="0",
                    on_change=self._on_change_handler,
                    text_align=ft.TextAlign.CENTER,
                    dense=True,
                    width=48,
                    height=36,
                    content_padding=4,
                    text_style=ft.TextStyle(
                        color=GUIColors.DIAGRAM_BLACK, size=font_size_default
                    ),
                )

    def _build_bp_table(
        self,
        font_size_default: int,
        font_size_micro: int,
        font_size_table_header: int,
    ) -> ft.Column:
        """Build and return the Base Pair Scores styled table container."""
        from amplifyp.dna import Nucleotides

        bp_columns = [
            ft.DataColumn(
                ft.Stack(
                    [
                        ft.canvas.Canvas(
                            [
                                ft.canvas.Line(
                                    0,
                                    8,
                                    70,
                                    34,
                                    paint=ft.Paint(
                                        color=GUIColors.DIVIDER_GREY,
                                        stroke_width=1,
                                    ),
                                )
                            ],
                            width=70,
                            height=36,
                        ),
                        ft.Container(
                            content=ft.Text(
                                "Template",
                                weight=ft.FontWeight.BOLD,
                                size=font_size_micro,
                            ),
                            alignment=ft.Alignment(1, -1),
                            padding=ft.Padding(0, 2, 5, 0),
                            width=70,
                            height=36,
                        ),
                        ft.Container(
                            content=ft.Text(
                                "Primer",
                                weight=ft.FontWeight.BOLD,
                                size=font_size_micro,
                            ),
                            alignment=ft.Alignment(-1, 1),
                            padding=ft.Padding(5, 0, 0, 2),
                            width=70,
                            height=36,
                        ),
                    ],
                    width=70,
                    height=36,
                )
            )
        ]

        for c_char in Nucleotides.TEMPLATE:
            if c_char == Nucleotides.GAP:
                continue
            bp_columns.append(
                ft.DataColumn(
                    ft.Container(
                        content=ft.Text(
                            c_char,
                            weight=ft.FontWeight.BOLD,
                            size=font_size_table_header,
                        ),
                        width=48,
                        alignment=ft.Alignment(0, 0),
                    )
                )
            )

        bp_rows = []
        for r_char in Nucleotides.PRIMER:
            cells = [
                ft.DataCell(
                    ft.Container(
                        content=ft.Text(
                            r_char,
                            weight=ft.FontWeight.BOLD,
                            size=font_size_table_header,
                        ),
                        width=70,
                        alignment=ft.Alignment(0, 0),
                    )
                )
            ]
            for c_char in Nucleotides.TEMPLATE:
                if c_char == Nucleotides.GAP:
                    continue
                key = f"bp_score_{r_char}_{c_char}"

                # Style TextField to match the borderless style in primer table
                self.settings_map[key].border = ft.InputBorder.NONE
                self.settings_map[key].width = 48
                self.settings_map[key].height = 30
                self.settings_map[key].content_padding = 0

                cells.append(
                    ft.DataCell(
                        ft.Container(
                            content=self.settings_map[key],
                            width=48,
                            alignment=ft.Alignment(0, 0),
                        )
                    )
                )

            bp_rows.append(ft.DataRow(cells=cells))

        bp_table = ft.DataTable(
            border=ft.Border.all(1, GUIColors.TRANSPARENT),
            vertical_lines=ft.BorderSide(1, GUIColors.DIVIDER_GREY),
            horizontal_lines=ft.BorderSide(1, GUIColors.DIVIDER_GREY),
            column_spacing=10,
            horizontal_margin=10,
            heading_row_color=GUIColors.INFO_HEADER_BG,
            heading_row_height=40,
            data_row_min_height=36,
            data_row_max_height=36,
            columns=bp_columns,
            rows=bp_rows,
        )

        return ft.Column(
            [
                ft.Text(
                    "Base Pair Weights",
                    weight=ft.FontWeight.BOLD,
                    size=font_size_default,
                ),
                ft.Row(
                    [
                        ft.Container(
                            content=bp_table,
                            border=ft.Border.all(1, GUIColors.OUTLINE),
                            border_radius=5,
                            padding=0,
                        )
                    ],
                    scroll=ft.ScrollMode.ADAPTIVE,
                    alignment=ft.MainAxisAlignment.CENTER,
                ),
            ],
            spacing=10,
            horizontal_alignment=ft.CrossAxisAlignment.CENTER,
        )

    def _build_pd_table(
        self,
        font_size_default: int,
        font_size_micro: int,
        font_size_table_header: int,
    ) -> ft.Column:
        """Build and return the Primer Dimer Scores styled table container."""
        from amplifyp.dna import Nucleotides

        pd_columns = [
            ft.DataColumn(
                ft.Stack(
                    [
                        ft.canvas.Canvas(
                            [
                                ft.canvas.Line(
                                    0,
                                    8,
                                    70,
                                    34,
                                    paint=ft.Paint(
                                        color=GUIColors.DIVIDER_GREY,
                                        stroke_width=1,
                                    ),
                                )
                            ],
                            width=70,
                            height=36,
                        ),
                        ft.Container(
                            content=ft.Text(
                                "Primer",
                                weight=ft.FontWeight.BOLD,
                                size=font_size_micro,
                            ),
                            alignment=ft.Alignment(1, -1),
                            padding=ft.Padding(0, 2, 5, 0),
                            width=70,
                            height=36,
                        ),
                        ft.Container(
                            content=ft.Text(
                                "Primer",
                                weight=ft.FontWeight.BOLD,
                                size=font_size_micro,
                            ),
                            alignment=ft.Alignment(-1, 1),
                            padding=ft.Padding(5, 0, 0, 2),
                            width=70,
                            height=36,
                        ),
                    ],
                    width=70,
                    height=36,
                )
            )
        ]

        for c_char in Nucleotides.PRIMER:
            pd_columns.append(
                ft.DataColumn(
                    ft.Container(
                        content=ft.Text(
                            c_char,
                            weight=ft.FontWeight.BOLD,
                            size=font_size_table_header,
                        ),
                        width=48,
                        alignment=ft.Alignment(0, 0),
                    )
                )
            )

        pd_rows = []
        for r_char in Nucleotides.PRIMER:
            cells = [
                ft.DataCell(
                    ft.Container(
                        content=ft.Text(
                            r_char,
                            weight=ft.FontWeight.BOLD,
                            size=font_size_table_header,
                        ),
                        width=70,
                        alignment=ft.Alignment(0, 0),
                    )
                )
            ]
            for c_char in Nucleotides.PRIMER:
                key = f"pd_score_{r_char}_{c_char}"

                # Style TextField to match the borderless style in primer table
                self.settings_map[key].border = ft.InputBorder.NONE
                self.settings_map[key].width = 48
                self.settings_map[key].height = 30
                self.settings_map[key].content_padding = 0

                cells.append(
                    ft.DataCell(
                        ft.Container(
                            content=self.settings_map[key],
                            width=48,
                            alignment=ft.Alignment(0, 0),
                        )
                    )
                )

            pd_rows.append(ft.DataRow(cells=cells))

        pd_table = ft.DataTable(
            border=ft.Border.all(1, GUIColors.TRANSPARENT),
            vertical_lines=ft.BorderSide(1, GUIColors.DIVIDER_GREY),
            horizontal_lines=ft.BorderSide(1, GUIColors.DIVIDER_GREY),
            column_spacing=10,
            horizontal_margin=10,
            heading_row_color=GUIColors.INFO_HEADER_BG,
            heading_row_height=40,
            data_row_min_height=36,
            data_row_max_height=36,
            columns=pd_columns,
            rows=pd_rows,
        )

        return ft.Column(
            [
                ft.Text(
                    "Primer Dimer Weights",
                    weight=ft.FontWeight.BOLD,
                    size=font_size_default,
                ),
                ft.Row(
                    [
                        ft.Container(
                            content=pd_table,
                            border=ft.Border.all(1, GUIColors.OUTLINE),
                            border_radius=5,
                            padding=0,
                        )
                    ],
                    scroll=ft.ScrollMode.ADAPTIVE,
                    alignment=ft.MainAxisAlignment.CENTER,
                ),
            ],
            spacing=10,
            horizontal_alignment=ft.CrossAxisAlignment.CENTER,
        )

    def _build_layout(
        self,
        bp_table_container: ft.Column,
        pd_table_container: ft.Column,
        header_size: int,
    ) -> None:
        """Construct self.controls layout list with ExpansionTiles."""
        self.controls = [
            self._build_origin_replication_tile(
                bp_table_container, header_size
            ),
            self._build_tm_settings_tile(header_size),
            self._build_primer_dimer_tile(pd_table_container, header_size),
            self._build_appearance_tile(header_size),
            ft.Divider(),
            self._build_action_buttons(),
        ]

        # Sync initial UI state
        self.update_ui()

    def _build_origin_replication_tile(
        self, bp_table_container: ft.Column, header_size: int
    ) -> ft.ExpansionTile:
        """Build the ExpansionTile for Origin of Replication settings."""
        return ft.ExpansionTile(
            title=ft.Text(
                "Origin of Replication Settings",
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Row(
                        [
                            bp_table_container,
                            ft.VerticalDivider(),
                            ft.Container(
                                content=ft.Column(
                                    [
                                        ft.Text(
                                            "Parameters",
                                            weight=ft.FontWeight.BOLD,
                                            size=self.settings.get(
                                                "font_size_default", 14
                                            ),
                                        ),
                                        self.set_primability_cutoff,
                                        self.set_stability_cutoff,
                                        self.set_amp4_compat,
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=220,
                            ),
                        ],
                        vertical_alignment=ft.CrossAxisAlignment.START,
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )

    def _build_tm_settings_tile(self, header_size: int) -> ft.ExpansionTile:
        """Build the ExpansionTile for Primer Melting Temperature settings."""
        return ft.ExpansionTile(
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

    def _build_primer_dimer_tile(
        self, pd_table_container: ft.Column, header_size: int
    ) -> ft.ExpansionTile:
        """Build the ExpansionTile for Primer Dimer settings."""
        return ft.ExpansionTile(
            title=ft.Text(
                "Primer Dimer Settings",
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Row(
                        [
                            pd_table_container,
                            ft.VerticalDivider(),
                            ft.Container(
                                content=ft.Column(
                                    [
                                        ft.Text(
                                            "Parameters",
                                            weight=ft.FontWeight.BOLD,
                                            size=self.settings.get(
                                                "font_size_default", 14
                                            ),
                                        ),
                                        self.set_pd_min_overlap,
                                        self.set_pd_threshold,
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=220,
                            ),
                        ],
                        vertical_alignment=ft.CrossAxisAlignment.START,
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )

    def _build_appearance_tile(self, header_size: int) -> ft.ExpansionTile:
        """Build the ExpansionTile for Appearance settings."""
        return ft.ExpansionTile(
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

    def _build_action_buttons(self) -> ft.Row:
        """Build the Action buttons Row (Apply & Reset)."""
        return ft.Row(
            [
                ft.FilledButton(
                    "Apply",
                    icon=ft.Icons.DONE,
                    on_click=self._on_apply_handler,
                ),
                ft.OutlinedButton(
                    "Reset to Default",
                    icon=ft.Icons.RESTORE,
                    on_click=self._on_reset_handler,
                ),
            ],
            alignment=ft.MainAxisAlignment.END,
            spacing=10,
        )

    def sync_to_state(self) -> None:
        """Sync current settings UI controls back to the central state."""
        for k, v in self.settings_map.items():
            self.settings[k] = v.value

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central settings."""
        for k, v in self.settings.items():
            if k in self.settings_map:
                if isinstance(self.settings_map[k], ft.Checkbox):
                    self.settings_map[k].value = bool(v)
                else:
                    self.settings_map[k].value = str(v)

    def _on_change_handler(self, e: ft.ControlEvent) -> None:
        """Handle change in settings fields."""
        self.sync_to_state()
        if self.on_change:
            self.on_change(e)

    def _on_apply_handler(self, e: ft.ControlEvent) -> None:
        """Handle apply button click."""
        self.sync_to_state()
        if self.on_apply:
            self.on_apply(e)

    def _on_reset_handler(self, e: ft.ControlEvent) -> None:
        """Handle reset to default button click."""
        from amplifyp.dna import Nucleotides
        from amplifyp.settings import (
            DEFAULT_BASE_PAIR_WEIGHTS,
            DEFAULT_PRIMABILITY_CUTOFF,
            DEFAULT_PRIMER_DIMER_OVERLAP,
            DEFAULT_PRIMER_DIMER_THRESHOLD,
            DEFAULT_PRIMER_DIMER_WEIGHTS,
            DEFAULT_STABILITY_CUTOFF,
            GLOBAL_TM_SETTINGS,
        )

        reset_dict = {
            "primability_cutoff": str(DEFAULT_PRIMABILITY_CUTOFF),
            "stability_cutoff": str(DEFAULT_STABILITY_CUTOFF),
            "amp4_compat": False,
            "tm_method": "SantaLucia 1998 / Owczarzy 2008 (Default)",
            "tm_dna_conc": str(GLOBAL_TM_SETTINGS.dna_conc),
            "tm_dnap_conc": str(GLOBAL_TM_SETTINGS.dnap_conc),
            "tm_mono_salt": str(GLOBAL_TM_SETTINGS.monovalent_salt_conc),
            "tm_div_salt": str(GLOBAL_TM_SETTINGS.divalent_salt_conc),
            "tm_dNTP_conc": str(GLOBAL_TM_SETTINGS.dnTP_conc),
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

        for r_char in Nucleotides.PRIMER:
            for c_char in Nucleotides.PRIMER:
                key = f"pd_score_{r_char}_{c_char}"
                reset_dict[key] = str(
                    DEFAULT_PRIMER_DIMER_WEIGHTS[r_char, c_char]
                )

        self.settings.from_dict(reset_dict)
        self.update_ui()
        self.app_page.update()
        if self.on_reset:
            self.on_reset(e)

    def get_replication_settings(self) -> ReplicationSettings:
        """Get the current settings as a ReplicationSettings object."""
        self.sync_to_state()
        return self.settings.get_replication_settings()

    def get_state(self) -> dict[str, Any]:
        """Get the current settings for serialization."""
        self.sync_to_state()
        return self.settings.to_dict()

    def set_state(self, state: dict[str, Any]) -> None:
        """Set the current settings from deserialized data."""
        self.settings.from_dict(state)
        self.update_ui()
        self.app_page.update()
