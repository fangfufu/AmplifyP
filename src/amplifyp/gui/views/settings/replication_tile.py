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

"""ReplicationTile expansion tile component for settings view."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.dna import Nucleotides
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.gui_helpers import (
    BorderedCheckbox,
    initialise_score_fields,
)
from amplifyp.gui.views.settings.lengthwise_weight_table import (
    LengthwiseWeightTable,
)
from amplifyp.gui.views.settings.score_table import ScoreTable
from amplifyp.settings import (
    DEFAULT_BASE_PAIR_WEIGHTS,
    DEFAULT_MATCH_WEIGHTS,
    DEFAULT_RUN_WEIGHTS,
)


class ReplicationTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for PCR settings."""

    def __init__(
        self,
        settings: GUISettings,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.Event | None], None],
        header_size: int,
        font_size_default: int,
        font_size_micro: int,
        font_size_table_header: int,
    ) -> None:
        """Initialise the ReplicationTile."""
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        self.set_primability_cutoff = ft.TextField(
            label="Primability Cutoff",
            value="0.8",
            expand=True,
            on_change=on_change_handler,
            border_color=GUIColours.OUTLINE,
        )
        self.set_stability_cutoff = ft.TextField(
            label="Stability Cutoff",
            value="0.4",
            expand=True,
            on_change=on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        self.set_amp4_compat = BorderedCheckbox(
            label="Amplify 4 Compatibility Mode",
            value=False,
            on_change=on_change_handler,
        )
        self.set_amp4_compat.expand = True
        self.set_improved_visualisation = BorderedCheckbox(
            label="Improved Replication Origin Visualisation",
            on_change=on_change_handler,
        )
        self.set_improved_visualisation.expand = True
        self.set_amplicon_ranking = ft.Dropdown(
            label="Amplicon vertical ranking",
            value=self.settings.get(
                "pcr_amplicon_ranking", "Position, then length"
            ),
            options=[
                ft.dropdown.Option("Position, then length"),
                ft.dropdown.Option("Quality score"),
                ft.dropdown.Option("Position, then quality"),
            ],
            expand=True,
            on_select=on_change_handler,
            text_size=font_size_default,
            border_color=GUIColours.OUTLINE,
        )
        self.set_primer_label_orientation = ft.Dropdown(
            label="Primer label orientation",
            value=self.settings.get(
                "pcr_primer_label_orientation", "Horizontal"
            ),
            options=[
                ft.dropdown.Option("Horizontal"),
                ft.dropdown.Option("45°"),
                ft.dropdown.Option("Vertical"),
            ],
            expand=True,
            on_select=on_change_handler,
            text_size=font_size_default,
            border_color=GUIColours.OUTLINE,
        )
        self.set_pcr_divider_behaviour = ft.Dropdown(
            label="Horizontal divider behaviour",
            value=self.settings.get("pcr_divider_behaviour", "autoresize"),
            options=[
                ft.dropdown.Option(
                    key="autoresize_1_4",
                    text="Auto-resize (max 1/4 of viewport)",
                ),
                ft.dropdown.Option(
                    key="autoresize",
                    text="Auto-resize (max 1/3 of viewport)",
                ),
                ft.dropdown.Option(
                    key="autoresize_1_2",
                    text="Auto-resize (max 1/2 of viewport)",
                ),
                ft.dropdown.Option(
                    key="autoresize_2_3",
                    text="Auto-resize (max 2/3 of viewport)",
                ),
                ft.dropdown.Option(
                    key="save_proportionality",
                    text="Save proportionality across session",
                ),
            ],
            expand=True,
            on_select=on_change_handler,
            text_size=font_size_default,
            border_color=GUIColours.OUTLINE,
        )

        settings_map["primability_cutoff"] = self.set_primability_cutoff
        settings_map["stability_cutoff"] = self.set_stability_cutoff
        settings_map["amp4_compat"] = self.set_amp4_compat
        settings_map["improved_visualisation"] = self.set_improved_visualisation
        settings_map["pcr_amplicon_ranking"] = self.set_amplicon_ranking
        settings_map["pcr_primer_label_orientation"] = (
            self.set_primer_label_orientation
        )
        settings_map["pcr_divider_behaviour"] = self.set_pcr_divider_behaviour

        col_headers = [c for c in Nucleotides.TEMPLATE if c != Nucleotides.GAP]

        initialise_score_fields(
            settings_map=self.settings_map,
            prefix="bp_score",
            row_headers=list(Nucleotides.PRIMER),
            col_headers=col_headers,
            on_change_handler=self.on_change_handler,
            font_size=font_size_default,
        )

        self.score_table = ScoreTable(
            label="Base Pair Weights",
            row_headers=list(Nucleotides.PRIMER),
            col_headers=col_headers,
            row_label="Primer",
            col_label="Template",
            prefix="bp_score",
            settings_map=self.settings_map,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
            scrollable=False,
            default_weights=DEFAULT_BASE_PAIR_WEIGHTS,
            on_change_handler=self.on_change_handler,
        )

        match_indices = sorted(
            int(k[len("match_weight_") :])
            for k in self.settings.keys()
            if k.startswith("match_weight_")
            and k[len("match_weight_") :].isdigit()
        )
        self.match_weight_table = LengthwiseWeightTable(
            label="Match Weights",
            indices=match_indices,
            prefix="match_weight",
            settings_map=self.settings_map,
            on_change_handler=self.on_change_handler,
            font_size_default=font_size_default,
            font_size_table_header=font_size_table_header,
            default_weights=DEFAULT_MATCH_WEIGHTS,
            row_header_label="3' Position",
        )

        run_indices = sorted(
            int(k[len("run_weight_") :])
            for k in self.settings.keys()
            if k.startswith("run_weight_") and k[len("run_weight_") :].isdigit()
        )
        self.run_weight_table = LengthwiseWeightTable(
            label="Run Weights",
            indices=run_indices,
            prefix="run_weight",
            settings_map=self.settings_map,
            on_change_handler=self.on_change_handler,
            font_size_default=font_size_default,
            font_size_table_header=font_size_table_header,
            default_weights=DEFAULT_RUN_WEIGHTS,
            row_header_label="Run Length",
        )

        super().__init__(
            title=ft.Text(
                "PCR",
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Column(
                        [
                            ft.Container(
                                content=ft.Column(
                                    [
                                        ft.Text(
                                            "Overview Map",
                                            weight=ft.FontWeight.BOLD,
                                            size=font_size_default + 2,
                                        ),
                                        ft.Row(
                                            [
                                                self.set_pcr_divider_behaviour,
                                                ft.Container(expand=True),
                                            ],
                                            spacing=15,
                                        ),
                                        ft.Divider(),
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=700,
                            ),
                            ft.Container(
                                content=ft.Text(
                                    "Replication Origin",
                                    weight=ft.FontWeight.BOLD,
                                    size=font_size_default + 2,
                                ),
                                width=700,
                            ),
                            ft.Row(
                                [
                                    self.score_table,
                                    self.match_weight_table,
                                    self.run_weight_table,
                                ],
                                wrap=True,
                                alignment=ft.MainAxisAlignment.CENTER,
                                vertical_alignment=ft.CrossAxisAlignment.START,
                                spacing=30,
                                run_spacing=20,
                            ),
                            ft.Container(height=10),
                            ft.Container(
                                content=ft.Column(
                                    [
                                        ft.Row(
                                            [
                                                self.set_primability_cutoff,
                                                self.set_stability_cutoff,
                                            ],
                                            spacing=15,
                                        ),
                                        ft.Row(
                                            [
                                                self.set_amp4_compat,
                                                self.set_improved_visualisation,
                                            ],
                                            spacing=15,
                                        ),
                                        ft.Row(
                                            [
                                                self.set_primer_label_orientation,
                                                ft.Container(expand=True),
                                            ],
                                            spacing=15,
                                        ),
                                        ft.Divider(),
                                        ft.Text(
                                            "Amplicon",
                                            weight=ft.FontWeight.BOLD,
                                            size=font_size_default + 2,
                                        ),
                                        ft.Row(
                                            [
                                                self.set_amplicon_ranking,
                                                ft.Container(expand=True),
                                            ],
                                            spacing=15,
                                        ),
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=700,
                            ),
                        ],
                        horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                        spacing=10,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )

    def update_ui(self) -> None:
        """Update Flet UI controls to match theme/settings."""
        self.match_weight_table.sync_from_settings(self.settings)
        self.run_weight_table.sync_from_settings(self.settings)
        self.score_table.table.heading_row_color = GUIColours.INFO_HEADER_BG
        self.match_weight_table.table.heading_row_color = (
            GUIColours.INFO_HEADER_BG
        )
        self.run_weight_table.table.heading_row_color = (
            GUIColours.INFO_HEADER_BG
        )
