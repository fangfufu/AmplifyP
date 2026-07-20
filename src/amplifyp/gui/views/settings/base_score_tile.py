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

"""Score table for Flet settings view."""

from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours


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
        width: int = 810,
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
            width: Total table container width in pixels. Defaults to 810.
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
                ft.canvas.Canvas(  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
                    [
                        ft.canvas.Line(  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
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

        columns = [ft.DataColumn(label=header_stack)]
        for i, c_char in enumerate(self.col_headers):
            is_last = i == len(self.col_headers) - 1
            col_width = 48 if is_last else 38
            columns.append(
                ft.DataColumn(
                    label=ft.Container(
                        content=ft.Text(
                            c_char,
                            weight=ft.FontWeight.BOLD,
                            size=self.font_size_table_header,
                        ),
                        width=col_width,
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
            for i, c_char in enumerate(self.col_headers):
                is_last = i == len(self.col_headers) - 1
                col_width = 48 if is_last else 38
                key = f"{self.prefix}_{r_char}_{c_char}"
                field = self.settings_map[key]

                # Style TextField to match the borderless style in primer table
                field.border = ft.InputBorder.NONE
                field.width = col_width
                field.height = 30
                field.content_padding = 0

                cells.append(
                    ft.DataCell(
                        ft.Container(
                            content=field,
                            width=col_width,
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
                            width=width,
                        )
                    ],
                    scroll=ft.ScrollMode.ALWAYS,
                    alignment=ft.MainAxisAlignment.CENTER,
                ),
            ],
            spacing=10,
            horizontal_alignment=ft.CrossAxisAlignment.CENTER,
        )
