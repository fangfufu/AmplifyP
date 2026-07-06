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

"""ListView for rendering, styling, and resizing primers."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import flet as ft

from .primer_row import PrimerRow
from .primer_validation import validate_primers

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .primer_input import PrimerInput


class PrimerList(ft.ListView):  # type: ignore[misc]
    """ListView for rendering, styling, and resizing primers."""

    def __init__(
        self,
        primer_input: PrimerInput,
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
            tooltip="Primer List",
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

        Reuses existing PrimerRow controls where possible, adjusting
        the list size incrementally. Updates highlights and header
        checkbox state.
        """
        font_family = self.primer_input.settings.get(
            "font_family", "Roboto Mono"
        )

        # Initialise with a single empty row if the list is completely empty
        if not self.primer_input.input_data.primers:
            self.primer_input.input_data.primers = [
                {"name": "", "seq": "", "active": False}
            ]

        ignore_inactive_name_dup = self.primer_input.settings.get(
            "ignore_inactive_name_dup_warn", True
        )
        ignore_inactive_seq_dup = self.primer_input.settings.get(
            "ignore_inactive_seq_dup_warn", True
        )

        self.primer_input.validation_errors = validate_primers(
            self.primer_input.input_data.primers,
            ignore_inactive_name_dup,
            ignore_inactive_seq_dup,
        )
        num_primers = len(self.primer_input.input_data.primers)

        # Clamp selection indices in case the list was truncated externally
        self.primer_input.selected_indices = {
            i for i in self.primer_input.selected_indices if i < num_primers
        }
        if (
            self.primer_input.focused_primer_index is not None
            and self.primer_input.focused_primer_index >= num_primers
        ):
            self.primer_input.focused_primer_index = None

        num_controls = len(self.controls)

        # 1. Adjust length of controls to match num_primers
        if num_controls > num_primers:
            self.controls[num_primers:] = []
        elif num_controls < num_primers:
            for idx in range(num_controls, num_primers):
                row = PrimerRow(
                    idx=idx,
                    name="",
                    seq="",
                    is_active=False,
                    is_dup=False,
                    name_error=None,
                    seq_error=None,
                    font_family=font_family,
                    name_column_width=self.primer_input.name_column_width,
                    settings=self.primer_input.settings,
                    on_change_handler=self.primer_input.on_change_handler,
                    handle_field_focus=self.primer_input.handle_field_focus,
                    handle_field_blur=self.primer_input.handle_field_blur,
                    handle_field_submit=self.primer_input.handle_field_submit,
                    on_row_click=self.primer_input.action_controller.handle_row_click,
                    on_row_double_click=self.primer_input.action_controller.handle_row_double_click,
                    on_divider_pan=self.primer_input.layout_manager.on_primer_divider_pan,
                    on_divider_pan_end=self.primer_input.layout_manager.on_primer_divider_pan_end,
                    is_focused=False,
                    is_last_row=(idx == num_primers - 1),
                    on_drag_start=self.primer_input.action_controller.handle_drag_start,
                    on_drag_update=self.primer_input.action_controller.handle_drag_update,
                    on_drag_end=self.primer_input.action_controller.handle_drag_end,
                )
                self.controls.append(row)

        # 2. Update controls in-place
        for idx, p in enumerate(self.primer_input.input_data.primers):
            name_val = p["name"]
            seq_val = p["seq"]
            is_active = p.get("active", True)
            error_message = self.primer_input.validation_errors[idx]

            row = self.controls[idx]
            if not isinstance(row, PrimerRow):
                continue

            row.update_index(
                idx,
                self.primer_input.action_controller.handle_row_click,
                self.primer_input.action_controller.handle_row_double_click,
            )

            if row.name_field.value != name_val:
                row.name_field.value = name_val
            if row.seq_field.value != seq_val:
                row.seq_field.value = seq_val
            if row.checkbox.value != is_active:
                row.checkbox.value = is_active

            row.set_error(error_message)
            row.update_tm(self.primer_input.settings)

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
                is_focused = c_idx in getattr(
                    self.primer_input, "selected_indices", set()
                )

                row.update_highlight_and_reorder(
                    is_focused=is_focused, is_dup=is_dup
                )
        try:
            if self.page:
                self.update()
        except RuntimeError:
            logger.debug("Page detached, skipping list update")
