# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""ScoreTable component for displaying base pair and primer dimer weights."""

from typing import Any

import flet as ft

from amplifyp.gui.settings import GUIColors


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
        """Initialize the ScoreTable."""
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
                        width=48,
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
                field.width = 48
                field.height = 30
                field.content_padding = 0

                cells.append(
                    ft.DataCell(
                        ft.Container(
                            content=field,
                            width=48,
                            alignment=ft.Alignment(0, 0),
                        )
                    )
                )

            rows.append(ft.DataRow(cells=cells))

        table = ft.DataTable(
            border=ft.Border.all(1, GUIColors.TRANSPARENT),
            vertical_lines=ft.BorderSide(1, GUIColors.DIVIDER_GREY),
            horizontal_lines=ft.BorderSide(1, GUIColors.DIVIDER_GREY),
            column_spacing=10,
            horizontal_margin=10,
            heading_row_color=GUIColors.INFO_HEADER_BG,
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
                            content=table,
                            border=ft.Border.all(1, GUIColors.OUTLINE),
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
