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

"""LengthwiseWeightTable helper component for settings view."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.settings import LengthWiseWeightTbl


class LengthwiseWeightTable(ft.Column):  # type: ignore[misc]
    """A Flet component displaying a vertical table of lengthwise weights."""

    def __init__(
        self,
        label: str,
        indices: list[int],
        prefix: str,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.Event | None], None],
        font_size_default: int,
        font_size_table_header: int,
        default_weights: LengthWiseWeightTbl,
    ) -> None:
        """Initialise the LengthwiseWeightTable."""
        self.label = label
        self.indices = sorted(indices)
        self.prefix = prefix
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler
        self.font_size_default = font_size_default
        self.font_size_table_header = font_size_table_header
        self.default_weights = default_weights

        for idx in self.indices:
            key = f"{self.prefix}_{idx}"
            if key not in self.settings_map:
                self.settings_map[key] = self._create_field("0")

        default_key = f"{self.prefix}_default"
        if default_key not in self.settings_map:
            self.settings_map[default_key] = self._create_field("0")

        self.add_button = ft.IconButton(
            icon=ft.Icons.ADD,
            tooltip="Add length entry",
            on_click=self._on_add_click,
            icon_size=16,
            padding=0,
            width=24,
            height=24,
        )
        self.remove_button = ft.IconButton(
            icon=ft.Icons.REMOVE,
            tooltip="Remove last length entry",
            on_click=self._on_remove_click,
            icon_size=16,
            padding=0,
            width=24,
            height=24,
        )
        self.reset_button = ft.TextButton(
            "Reset",
            icon=ft.Icons.RESTORE,
            on_click=self._on_reset_click,
        )

        columns = [
            ft.DataColumn(
                label=ft.Container(
                    content=ft.Text(
                        "Length",
                        weight=ft.FontWeight.BOLD,
                        size=self.font_size_table_header,
                    ),
                    width=60,
                    alignment=ft.Alignment(0, 0),
                )
            ),
            ft.DataColumn(
                label=ft.Container(
                    content=ft.Text(
                        "Weight",
                        weight=ft.FontWeight.BOLD,
                        size=self.font_size_table_header,
                    ),
                    width=65,
                    alignment=ft.Alignment(0, 0),
                )
            ),
        ]

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
            rows=[],
        )
        self._rebuild_rows()

        super().__init__(
            controls=[
                ft.Text(
                    self.label,
                    weight=ft.FontWeight.BOLD,
                    size=self.font_size_default,
                ),
                ft.Container(
                    content=self.table,
                    border=ft.Border.all(1, GUIColours.OUTLINE),
                    border_radius=5,
                    padding=0,
                ),
                ft.Row(
                    [
                        self.add_button,
                        self.remove_button,
                        self.reset_button,
                    ],
                    alignment=ft.MainAxisAlignment.CENTER,
                    spacing=4,
                    tight=True,
                ),
            ],
            spacing=8,
            horizontal_alignment=ft.CrossAxisAlignment.CENTER,
        )

    def _create_field(self, value: str) -> ft.TextField:
        """Create a borderless TextField matching ScoreTable styling."""
        return ft.TextField(
            value=value,
            width=65,
            height=30,
            content_padding=0,
            dense=True,
            border=ft.NoInputBorder(),
            text_align=ft.TextAlign.CENTER,
            text_vertical_align=ft.VerticalAlignment.CENTER,
            text_size=self.font_size_default,
            on_change=self.on_change_handler,
        )

    def _make_row(self, row_label: str, field: ft.TextField) -> ft.DataRow:
        """Create a DataRow with a centred label and weight field."""
        return ft.DataRow(
            cells=[
                ft.DataCell(
                    ft.Container(
                        content=ft.Text(
                            row_label,
                            weight=ft.FontWeight.BOLD,
                            size=self.font_size_table_header,
                        ),
                        width=60,
                        alignment=ft.Alignment(0, 0),
                    )
                ),
                ft.DataCell(
                    ft.Container(
                        content=field,
                        width=65,
                        alignment=ft.Alignment(0, 0),
                    )
                ),
            ]
        )

    def _rebuild_rows(self) -> None:
        """Rebuild the DataTable rows from current indices and Default."""
        rows = [
            self._make_row(str(idx), self.settings_map[f"{self.prefix}_{idx}"])
            for idx in self.indices
        ]
        rows.append(
            self._make_row(
                "Default", self.settings_map[f"{self.prefix}_default"]
            )
        )
        self.table.rows = rows

    def _on_add_click(self, _e: ft.Event | None) -> None:
        """Add the next length entry."""
        next_idx = (self.indices[-1] + 1) if self.indices else 0
        default_val = str(
            self.settings_map[f"{self.prefix}_default"].value or "0"
        )
        self.settings_map[f"{self.prefix}_{next_idx}"] = self._create_field(
            default_val
        )
        self.indices.append(next_idx)
        self._rebuild_rows()
        self.on_change_handler(None)
        if self.page:
            self.update()

    def _on_remove_click(self, _e: ft.Event | None) -> None:
        """Remove the last length entry, always keeping Default."""
        if not self.indices:
            return
        last_idx = self.indices.pop()
        self.settings_map.pop(f"{self.prefix}_{last_idx}", None)
        self._rebuild_rows()
        self.on_change_handler(None)
        if self.page:
            self.update()

    def _on_reset_click(self, _e: ft.Event | None) -> None:
        """Reset this table's entries to their default values."""
        for idx in self.indices:
            self.settings_map.pop(f"{self.prefix}_{idx}", None)

        self.indices = sorted(self.default_weights.overrides.keys())
        for idx in self.indices:
            val = self.default_weights.overrides[idx]
            self.settings_map[f"{self.prefix}_{idx}"] = self._create_field(
                str(val)
            )

        self.settings_map[f"{self.prefix}_default"].value = str(
            self.default_weights.default_weight
        )
        self._rebuild_rows()
        self.on_change_handler(None)
        if self.page:
            self.update()
