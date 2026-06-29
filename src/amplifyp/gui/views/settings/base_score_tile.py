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

"""Base score settings tile and score table for Flet settings view."""

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings


class ScoreTable(ft.Column):  # type: ignore[misc]
    """A Flet component displaying a grid/table of weight TextFields."""

    def __init__(
        self,
        label: str,
        row_headers: list[str],
        col_headers: list[str],
        row_label: str,
        col_label: str,
        prefix: str,
        settings_map: dict[str, Any],
        font_size_default: int,
        font_size_micro: int,
        font_size_table_header: int,
    ) -> None:
        """Initialise the ScoreTable.

        Args:
            label: The title text displayed above the table.
            row_headers: List of row header labels.
            col_headers: List of column header labels.
            row_label: Label for the row header column.
            col_label: Label for the column header row.
            prefix: Key prefix used to look up TextField components in
                settings_map (e.g. "pd_score", "bp_score").
            settings_map: Dictionary mapping setting keys to UI components.
            font_size_default: Default font size for the label.
            font_size_micro: Font size for the diagonal header labels.
            font_size_table_header: Font size for row and column headers.
        """
        self.label = label
        self.row_headers = row_headers
        self.col_headers = col_headers
        self.row_label = row_label
        self.col_label = col_label
        self.prefix = prefix
        self.settings_map = settings_map
        self.font_size_default = font_size_default
        self.font_size_micro = font_size_micro
        self.font_size_table_header = font_size_table_header

        # Diagonal line canvas header
        header_stack = ft.Stack(
            [
                ft.canvas.Canvas(
                    [
                        ft.canvas.Line(
                            0,
                            8,
                            70,
                            34,
                            paint=ft.Paint(
                                color=GUIColours.DIVIDER_GREY,
                                stroke_width=1,
                            ),
                        )
                    ],
                    width=70,
                    height=36,
                ),
                ft.Container(
                    content=ft.Text(
                        self.col_label,
                        weight=ft.FontWeight.BOLD,
                        size=self.font_size_micro,
                    ),
                    alignment=ft.Alignment(1, -1),
                    padding=ft.Padding(0, 2, 5, 0),
                    width=70,
                    height=36,
                ),
                ft.Container(
                    content=ft.Text(
                        self.row_label,
                        weight=ft.FontWeight.BOLD,
                        size=self.font_size_micro,
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

        columns = [ft.DataColumn(header_stack)]
        for c_char in self.col_headers:
            columns.append(
                ft.DataColumn(
                    ft.Container(
                        content=ft.Text(
                            c_char,
                            weight=ft.FontWeight.BOLD,
                            size=self.font_size_table_header,
                        ),
                        width=38,
                        alignment=ft.Alignment(0, 0),
                    )
                )
            )

        rows = []
        for r_char in self.row_headers:
            cells = [
                ft.DataCell(
                    ft.Container(
                        content=ft.Text(
                            r_char,
                            weight=ft.FontWeight.BOLD,
                            size=self.font_size_table_header,
                        ),
                        width=70,
                        alignment=ft.Alignment(0, 0),
                    )
                )
            ]
            for c_char in self.col_headers:
                key = f"{self.prefix}_{r_char}_{c_char}"
                field = self.settings_map[key]

                # Style TextField to match the borderless style in primer table
                field.border = ft.InputBorder.NONE
                field.width = 38
                field.height = 30
                field.content_padding = 0

                cells.append(
                    ft.DataCell(
                        ft.Container(
                            content=field,
                            width=38,
                            alignment=ft.Alignment(0, 0),
                        )
                    )
                )

            rows.append(ft.DataRow(cells=cells))

        self.table = ft.DataTable(
            border=ft.Border.all(1, GUIColours.TRANSPARENT),
            vertical_lines=ft.BorderSide(1, GUIColours.DIVIDER_GREY),
            horizontal_lines=ft.BorderSide(1, GUIColours.DIVIDER_GREY),
            column_spacing=10,
            horizontal_margin=10,
            heading_row_color=GUIColours.INFO_HEADER_BG,
            heading_row_height=40,
            data_row_min_height=36,
            data_row_max_height=36,
            columns=columns,
            rows=rows,
        )

        super().__init__(
            [
                ft.Text(
                    self.label,
                    weight=ft.FontWeight.BOLD,
                    size=self.font_size_default,
                ),
                ft.Row(
                    [
                        ft.Container(
                            content=self.table,
                            border=ft.Border.all(1, GUIColours.OUTLINE),
                            border_radius=5,
                            padding=0,
                        )
                    ],
                    scroll=ft.ScrollMode.ALWAYS,
                    alignment=ft.MainAxisAlignment.CENTER,
                ),
            ],
            spacing=10,
            horizontal_alignment=ft.CrossAxisAlignment.CENTER,
        )


class BaseScoreTile(ft.ExpansionTile):  # type: ignore[misc]
    """Base settings tile with a ScoreTable and parameter controls."""

    def __init__(
        self,
        settings: GUISettings,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.Event | None], None],
        header_size: int,
        font_size_default: int,
        font_size_micro: int,
        font_size_table_header: int,
        title: str,
        score_table_label: str,
        score_table_prefix: str,
        row_headers: list[str],
        col_headers: list[str],
        row_label: str,
        col_label: str,
        parameter_controls: list[ft.Control],
    ) -> None:
        """Initialise the BaseScoreTile.

        Args:
            settings: The settings object.
            settings_map: A dictionary mapping setting keys to UI
                components for population and retrieval.
            on_change_handler: The handler to call when a setting changes.
            header_size: The size of the expansion tile header text.
            font_size_default: Default font size for text elements.
            font_size_micro: Micro font size for small labels.
            font_size_table_header: Font size for table header cells.
            title: The expansion tile title text.
            score_table_label: Label for the score table section.
            score_table_prefix: Key prefix for score table fields in
                settings_map.
            row_headers: List of row header labels for the score table.
            col_headers: List of column header labels for the score table.
            row_label: Label for the row header column.
            col_label: Label for the column header row.
            parameter_controls: List of additional parameter controls to
                display alongside the score table.
        """
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        # Map initialisation (must happen before building ScoreTable)
        from amplifyp.gui.util import initialise_score_fields

        initialise_score_fields(
            settings_map=self.settings_map,
            prefix=score_table_prefix,
            row_headers=row_headers,
            col_headers=col_headers,
            on_change_handler=self.on_change_handler,
            font_size=font_size_default,
        )

        self.score_table = ScoreTable(
            label=score_table_label,
            row_headers=row_headers,
            col_headers=col_headers,
            row_label=row_label,
            col_label=col_label,
            prefix=score_table_prefix,
            settings_map=self.settings_map,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
        )

        super().__init__(
            title=ft.Text(
                title,
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
                                        *parameter_controls,
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=350,
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
