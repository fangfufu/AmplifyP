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

"""ReplicationTile component for Flet settings view."""

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.views.settings.base_score_tile import ScoreTable


class ReplicationTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Origin of Replication settings."""

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
        """Initialise the ReplicationTile.

        Args:
            settings (Any): The settings object.
            settings_map (dict[str, Any]): A dictionary mapping setting keys
                to UI components.
            on_change_handler (Any): The handler to call when a setting changes.
            header_size (int): The size of the header text.
            font_size_default (int): Default font size for text.
            font_size_micro (int): Micro font size for small text.
            font_size_table_header (int): Font size for table headers.
        """
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        from amplifyp.dna import Nucleotides
        from amplifyp.gui.utils.ui import (
            BorderedCheckbox,
            initialise_score_fields,
        )

        self.set_primability_cutoff = ft.TextField(
            label="Primability Cutoff",
            value="0.8",
            on_change=on_change_handler,
            border_color=GUIColours.OUTLINE,
        )
        self.set_stability_cutoff = ft.TextField(
            label="Stability Cutoff",
            value="0.4",
            on_change=on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        self.set_amp4_compat = BorderedCheckbox(
            label="Amplify4 Compatibility Mode",
            value=False,
            on_change=on_change_handler,
        )
        self.set_improved_visualisation = BorderedCheckbox(
            label="Improved Primer Binding Site Visualisation",
            on_change=on_change_handler,
        )

        settings_map["primability_cutoff"] = self.set_primability_cutoff
        settings_map["stability_cutoff"] = self.set_stability_cutoff
        settings_map["amp4_compat"] = self.set_amp4_compat
        settings_map["improved_visualisation"] = self.set_improved_visualisation

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
        )

        super().__init__(
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
                            self.score_table,
                            ft.Container(height=10),
                            ft.Container(
                                content=ft.Column(
                                    [
                                        ft.Divider(),
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
                                        self.set_improved_visualisation,
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=500,
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
        self.score_table.table.heading_row_color = GUIColours.INFO_HEADER_BG
