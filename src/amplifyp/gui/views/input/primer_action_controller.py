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

"""Action mutation controller for DNA primers input."""

import inspect
from typing import TYPE_CHECKING

import flet as ft

if TYPE_CHECKING:
    from .primer_input import PrimerInput


class PrimerActionController:
    """Handles mutation operations on the primers list.

    Such as adding, moving, and deleting.
    """

    def __init__(self, owner: "PrimerInput") -> None:
        """Initialise the PrimerActionController."""
        self.owner = owner

    def handle_row_click(self, idx: int, name_edit: ft.TextField) -> None:
        """Handle clicking on the row container.

        Selects the row and focuses the name field.
        """
        self.owner.focused_primer_index = idx
        self.owner._update_row_highlights()
        self.owner._update_primer_info_panel()
        res = name_edit.focus()
        if inspect.iscoroutine(res):
            if self.owner.app_page:

                async def do_focus() -> None:
                    await res

                self.owner.app_page.run_task(do_focus)

    def on_add_primer_row(self, idx: int) -> None:
        """Add a new empty primer row immediately below the row at idx."""
        self.owner.sync_to_state(rebuild_if_needed=False)
        self.owner.input_data.primers.insert(
            idx + 1, {"name": "", "seq": "", "active": False}
        )
        self.owner.update_ui()
        if self.owner.on_change_handler:
            self.owner.on_change_handler(None)

    def move_primer(self, idx: int, direction: int) -> None:
        """Move primer at idx up (direction=-1) or down (direction=1)."""
        self.owner.sync_to_state(rebuild_if_needed=False)
        primers = self.owner.input_data.primers
        target_idx = idx + direction

        # Swap only if both indices are valid filled primers.
        if 0 <= idx < len(primers) and 0 <= target_idx < len(primers):
            primers[idx], primers[target_idx] = (
                primers[target_idx],
                primers[idx],
            )
            if self.owner.focused_primer_index == idx:
                self.owner.focused_primer_index = target_idx
            elif self.owner.focused_primer_index == target_idx:
                self.owner.focused_primer_index = idx

            self.owner.update_ui()
            if self.owner.on_change_handler:
                self.owner.on_change_handler(None)

    def delete_primers(self, indices_to_delete: set[int]) -> None:
        """Delete primers at indices and re-index controls in-place."""
        if not indices_to_delete:
            return

        self.owner.sync_to_state(rebuild_if_needed=False)
        primers = self.owner.input_data.primers

        # Keep only indices NOT in the deleted set
        new_primers = [
            p for i, p in enumerate(primers) if i not in indices_to_delete
        ]
        if not new_primers:
            new_primers = [{"name": "", "seq": "", "active": False}]
        self.owner.input_data.primers = new_primers

        # Adjust focus index
        if self.owner.focused_primer_index in indices_to_delete:
            self.owner.focused_primer_index = None
        elif self.owner.focused_primer_index is not None:
            # Shift focus index down for each deleted index that was before it
            shift = sum(
                1
                for i in indices_to_delete
                if i < self.owner.focused_primer_index
            )
            self.owner.focused_primer_index -= shift

        from .primer_row import PrimerRow

        # Filter the controls list
        remaining_controls = []
        for row in self.owner.primers_list.controls:
            if isinstance(row, PrimerRow):
                if row.data in indices_to_delete:
                    continue
                remaining_controls.append(row)

        if not remaining_controls:
            self.owner.primers_list.controls.clear()
            self.owner.update_ui()
            return

        self.owner.primers_list.controls = remaining_controls
        num_remaining = len(remaining_controls)
        # Find the earliest deleted index to start re-indexing from
        start_reindex = min(indices_to_delete)

        for new_idx in range(start_reindex, num_remaining):
            row = remaining_controls[new_idx]
            if isinstance(row, PrimerRow):
                is_last_row = new_idx == num_remaining - 1
                row.update_index(
                    new_idx=new_idx,
                    is_last_row=is_last_row,
                    on_move_primer=self.move_primer,
                    on_delete_primer=lambda idx: self.delete_primers({idx}),
                    on_add_primer=self.on_add_primer_row,
                    on_row_click=self.handle_row_click,
                )

        self.owner._update_row_highlights()
        self.owner._update_primer_info_panel()
        self.owner._update_header_checkbox_state()
        self.owner._update_delete_button_disabled_state()
        if self.owner.app_page:
            self.owner.app_page.update()

        if self.owner.on_change_handler:
            self.owner.on_change_handler(None)
