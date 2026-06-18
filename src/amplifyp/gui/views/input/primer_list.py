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

"""ListView for rendering, styling, and resizing primers."""

from typing import Any

import flet as ft

from .primer_row import PrimerRow
from .primer_validation import validate_primers


class PrimerList(ft.ListView):  # type: ignore[misc]
    """ListView for rendering, styling, and resizing primers."""

    def __init__(
        self,
        primer_input: Any,
    ) -> None:
        """Initialise the PrimerList.

        Args:
            primer_input: The parent PrimerInput component that owns this
                list.
        """
        super().__init__(
            expand=True,
            spacing=0,
            padding=0,
            scroll=ft.ScrollMode.ALWAYS,
            on_scroll=self._on_scroll,
        )
        self.primer_input = primer_input
        self.scroll_pixels = 0.0
        self.viewport_dimension = 600.0

    def _on_scroll(self, e: ft.OnScrollEvent) -> None:
        """Update the cached scroll position and viewport dimension.

        Args:
            e: The Flet scroll event containing pixel position.
        """
        self.scroll_pixels = e.pixels

    def update_list_ui(self) -> None:
        """Update Flet UI controls to match the central state.

        Clears existing controls, creates a new PrimerRow for each
        primer in the input data, validates all primers, and updates
        highlights and header checkbox state.
        """
        font_family = self.primer_input.settings.get(
            "font_family", "Roboto Mono"
        )
        self.controls.clear()

        # Initialise with a single empty row if the list is completely empty
        if not self.primer_input.input_data.primers:
            self.primer_input.input_data.primers = [
                {"name": "", "seq": "", "active": False}
            ]

        dup_indices = self.primer_input._get_duplicate_indices()

        self.primer_input.validation_errors = validate_primers(
            self.primer_input.input_data.primers
        )
        num_primers = len(self.primer_input.input_data.primers)
        for idx, p in enumerate(self.primer_input.input_data.primers):
            name_val = p["name"]
            seq_val = p["seq"]
            is_active = p.get("active", True)
            error_message = self.primer_input.validation_errors[idx]
            name_err = (
                error_message.get("name")
                if isinstance(error_message, dict)
                else (
                    error_message
                    if error_message == "Duplicate primer name"
                    else None
                )
            )
            seq_err = (
                error_message.get("seq")
                if isinstance(error_message, dict)
                else (
                    None
                    if error_message == "Duplicate primer name"
                    else error_message
                )
            )

            is_dup = idx in dup_indices
            is_focused = idx == self.primer_input.focused_primer_index
            is_last_row = idx == num_primers - 1

            row = PrimerRow(
                idx=idx,
                name=name_val,
                seq=seq_val,
                is_active=is_active,
                is_dup=is_dup,
                name_error=name_err,
                seq_error=seq_err,
                font_family=font_family,
                name_column_width=self.primer_input.name_column_width,
                on_change_handler=self.primer_input.on_change_handler,
                handle_field_focus=self.primer_input.handle_field_focus,
                handle_field_blur=self.primer_input.handle_field_blur,
                handle_field_submit=self.primer_input.handle_field_submit,
                on_row_click=self.primer_input.action_controller.handle_row_click,
                on_move_primer=self.primer_input.action_controller.move_primer,
                on_delete_primer=lambda idx: (
                    self.primer_input.action_controller.delete_primers({idx})
                ),
                on_add_primer=self.primer_input.action_controller.on_add_primer_row,
                on_divider_pan=self.primer_input.layout_manager.on_primer_divider_pan,
                on_divider_pan_end=self.primer_input.layout_manager.on_primer_divider_pan_end,
                is_focused=is_focused,
                is_last_row=is_last_row,
            )
            self.controls.append(row)

        self.update_row_highlights()
        self.primer_input._update_primer_info_panel()
        self.primer_input._update_header_checkbox_state()

    def update_row_highlights(self) -> None:
        """Update background colours of all row containers.

        Highlights rows based on selection (focused primer) and
        duplicates (by name or sequence).
        """
        dup_indices = self.primer_input._get_duplicate_indices()

        for row in self.controls:
            if isinstance(row, PrimerRow) and row.data is not None:
                c_idx = row.data
                is_dup = c_idx in dup_indices
                is_focused = c_idx == self.primer_input.focused_primer_index
                row.update_highlight_and_reorder(
                    is_focused=is_focused, is_dup=is_dup
                )
        try:
            if self.page:
                self.update()
        except RuntimeError:
            pass
