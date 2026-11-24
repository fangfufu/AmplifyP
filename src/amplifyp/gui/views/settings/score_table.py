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

"""ScoreTable helper component for settings view."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.views.settings.weight_table_base import WeightTableBase
from amplifyp.settings import BasePairWeightsTbl


class ScoreTable(WeightTableBase):
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
        width: int | None = None,
        scrollable: bool = True,
        default_weights: BasePairWeightsTbl | None = None,
        on_change_handler: Callable[[ft.Event | None], None] | None = None,
    ) -> None:
        """Initialise the ScoreTable."""
        super().__init__(
            label=label,
            prefix=prefix,
            settings_map=settings_map,
            font_size_default=font_size_default,
            font_size_table_header=font_size_table_header,
            on_change_handler=on_change_handler,
        )
        self.row_headers = row_headers
        self.col_headers = col_headers
        self.row_label = row_label
        self.col_label = col_label
        self.font_size_micro = font_size_micro
        self.default_weights = default_weights

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

        col_width = 30
        columns = [ft.DataColumn(label=header_stack)]
        for c_char in self.col_headers:
            columns.append(
                ft.DataColumn(
                    label=self._build_header_label(c_char, width=col_width)
                )
            )

        rows = []
        for r_char in self.row_headers:
            cells = [ft.DataCell(self._build_header_label(r_char, width=70))]
            for c_char in self.col_headers:
                key = f"{self.prefix}_{r_char}_{c_char}"
                field = self._style_field(
                    self.settings_map[key], width=col_width
                )
                cells.append(
                    ft.DataCell(self._wrap_centred(field, width=col_width))
                )

            rows.append(ft.DataRow(cells=cells))

        self.table = self._create_data_table(
            columns=columns,
            rows=rows,
            column_spacing=4,
            horizontal_margin=4,
        )

        table_container = self._create_table_container(self.table, width=width)
        table_control: ft.Control = (
            ft.Row(
                [table_container],
                scroll=ft.ScrollMode.ALWAYS,
                alignment=ft.MainAxisAlignment.CENTER,
                tight=True,
            )
            if scrollable
            else table_container
        )

        self.reset_button: ft.TextButton | None = None
        action_controls: list[ft.Control] | None = None
        if self.default_weights is not None:
            self.reset_button = self._create_reset_button(self._on_reset_click)
            action_controls = [self.reset_button]

        self._init_column(
            table_control=table_control,
            action_controls=action_controls,
            spacing=10,
        )

    def _on_reset_click(self, _e: ft.Event | None) -> None:
        """Reset this table's entries to their default values."""
        if self.default_weights is None:
            return
        for r_char in self.row_headers:
            for c_char in self.col_headers:
                key = f"{self.prefix}_{r_char}_{c_char}"
                if key in self.settings_map:
                    self.settings_map[key].value = str(
                        self.default_weights[r_char, c_char]
                    )
        self._notify_changed()
