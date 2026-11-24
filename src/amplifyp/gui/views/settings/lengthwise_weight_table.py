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

from amplifyp.gui.settings import GUISettings
from amplifyp.gui.views.settings.weight_table_base import WeightTableBase
from amplifyp.settings import LengthWiseWeightTbl


class LengthwiseWeightTable(WeightTableBase):
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
        row_header_label: str = "Length",
    ) -> None:
        """Initialise the LengthwiseWeightTable."""
        super().__init__(
            label=label,
            prefix=prefix,
            settings_map=settings_map,
            font_size_default=font_size_default,
            font_size_table_header=font_size_table_header,
            on_change_handler=on_change_handler,
            spacing=8,
        )
        self.indices = sorted(indices)
        self.default_weights = default_weights
        self.row_header_label = row_header_label
        self._label_col_width = 92

        for idx in self.indices:
            key = f"{self.prefix}_{idx}"
            if key not in self.settings_map:
                self.settings_map[key] = self._create_field("0")

        default_key = f"{self.prefix}_default"
        if default_key not in self.settings_map:
            self.settings_map[default_key] = self._create_field("0")

        entry_noun = self.row_header_label.lower()
        self.add_button = ft.IconButton(
            icon=ft.Icons.ADD,
            tooltip=f"Add {entry_noun} entry",
            on_click=self._on_add_click,
            icon_size=16,
            padding=0,
            width=24,
            height=24,
        )
        self.remove_button = ft.IconButton(
            icon=ft.Icons.REMOVE,
            tooltip=f"Remove last {entry_noun} entry",
            on_click=self._on_remove_click,
            icon_size=16,
            padding=0,
            width=24,
            height=24,
        )
        self.reset_button = self._create_reset_button(self._on_reset_click)

        columns = [
            ft.DataColumn(
                label=self._build_header_label(
                    self.row_header_label, width=self._label_col_width
                )
            ),
            ft.DataColumn(label=self._build_header_label("Weight", width=65)),
        ]

        self.table = self._create_data_table(
            columns=columns,
            rows=[],
            column_spacing=10,
            horizontal_margin=10,
        )
        self._rebuild_rows()

        self._init_column(
            table_control=self._create_table_container(self.table),
            action_controls=[
                self.add_button,
                self.remove_button,
                self.reset_button,
            ],
            spacing=8,
        )

    def _create_field(self, value: str) -> ft.TextField:
        """Create a borderless TextField matching ScoreTable styling."""
        return self._create_cell_field(value, width=65)

    def _make_row(self, row_label: str, field: ft.TextField) -> ft.DataRow:
        """Create a DataRow with a centred label and weight field."""
        return ft.DataRow(
            cells=[
                ft.DataCell(
                    self._build_header_label(
                        row_label, width=self._label_col_width
                    )
                ),
                ft.DataCell(self._wrap_centred(field, width=65)),
            ]
        )

    def _rebuild_rows(self) -> None:
        """Rebuild the DataTable rows from current indices and Default."""
        rows = [
            self._make_row(
                str(idx + 1), self.settings_map[f"{self.prefix}_{idx}"]
            )
            for idx in self.indices
        ]
        rows.append(
            self._make_row(
                "Default", self.settings_map[f"{self.prefix}_default"]
            )
        )
        self.table.rows = rows

    def sync_from_settings(self, settings: GUISettings) -> None:
        """Sync rows and field values from GUISettings."""
        prefix_with_sep = f"{self.prefix}_"
        new_indices = sorted(
            int(k[len(prefix_with_sep) :])
            for k in settings.keys()
            if k.startswith(prefix_with_sep)
            and k[len(prefix_with_sep) :].isdigit()
        )
        if new_indices != self.indices:
            for idx in self.indices:
                self.settings_map.pop(f"{self.prefix}_{idx}", None)
            self.indices = new_indices
            for idx in self.indices:
                key = f"{self.prefix}_{idx}"
                self.settings_map[key] = self._create_field(
                    str(settings.get(key, self.default_weights.default_weight))
                )
            self._rebuild_rows()
        else:
            for idx in self.indices:
                key = f"{self.prefix}_{idx}"
                val = str(
                    settings.get(key, self.default_weights.default_weight)
                )
                if self.settings_map[key].value != val:
                    self.settings_map[key].value = val

        default_key = f"{self.prefix}_default"
        default_val = str(
            settings.get(default_key, self.default_weights.default_weight)
        )
        if self.settings_map[default_key].value != default_val:
            self.settings_map[default_key].value = default_val

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
        self._notify_changed()

    def _on_remove_click(self, _e: ft.Event | None) -> None:
        """Remove the last length entry, always keeping Default."""
        if not self.indices:
            return
        last_idx = self.indices.pop()
        self.settings_map.pop(f"{self.prefix}_{last_idx}", None)
        self._rebuild_rows()
        self._notify_changed()

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
        self._notify_changed()
