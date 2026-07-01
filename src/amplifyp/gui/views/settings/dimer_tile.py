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

"""DimerTile component for Flet settings view."""

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.views.settings.base_score_tile import ScoreTable


class DimerTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Dimer settings."""

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
        """Initialise the DimerTile.

        Args:
            settings: The settings object.
            settings_map: A dictionary mapping setting keys to UI
                components for population and retrieval.
            on_change_handler: The handler to call when a setting changes.
            header_size: The size of the expansion tile header text.
            font_size_default: Default font size for text elements.
            font_size_micro: Micro font size for small labels.
            font_size_table_header: Font size for table header cells.
        """
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        from amplifyp.dna import Nucleotides
        from amplifyp.gui.utils.ui import initialise_score_fields

        self.set_pd_min_overlap = ft.TextField(
            label="Min Overlap",
            value="3",
            on_change=on_change_handler,
            border_color=GUIColours.OUTLINE,
        )
        self.set_pd_threshold = ft.TextField(
            label="Threshold",
            value="60.0",
            on_change=on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        settings_map["pd_min_overlap"] = self.set_pd_min_overlap
        settings_map["pd_threshold"] = self.set_pd_threshold

        initialise_score_fields(
            settings_map=self.settings_map,
            prefix="pd_score",
            row_headers=list(Nucleotides.PRIMER),
            col_headers=list(Nucleotides.PRIMER),
            on_change_handler=self.on_change_handler,
            font_size=font_size_default,
        )

        self.score_table = ScoreTable(
            label="Primer Dimer Weights",
            row_headers=list(Nucleotides.PRIMER),
            col_headers=list(Nucleotides.PRIMER),
            row_label="Primer",
            col_label="Primer",
            prefix="pd_score",
            settings_map=self.settings_map,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
        )

        super().__init__(
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
                            self.score_table,
                            ft.Container(height=10),
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
