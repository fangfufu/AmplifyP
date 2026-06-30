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
        """Initialise the PrimerActionController.

        Args:
            owner: The parent PrimerInput component that owns this controller.
        """
        self.owner = owner
        self.current_drag_y = 0.0
        self.drag_start_index = 0

    def handle_row_click(self, idx: int, name_edit: ft.TextField) -> None:
        """Handle clicking on the row container.

        Selects the row, focuses the name field, updates highlights,
        and displays primer info.

        Args:
            idx: Zero-based index of the clicked primer row.
            name_edit: The name TextField control to focus.
        """
        if idx in self.owner.selected_indices:
            self.owner.selected_indices.remove(idx)
            if self.owner.focused_primer_index == idx:
                self.owner.focused_primer_index = None
        else:
            self.owner.selected_indices.add(idx)
            self.owner.focused_primer_index = idx

            res = name_edit.focus()
            if inspect.iscoroutine(res) and self.owner.app_page:

                async def do_focus() -> None:
                    await res

                self.owner.app_page.run_task(do_focus)

        self.owner._update_row_highlights()
        self.owner._update_primer_info_panel()
        self.owner._update_delete_button_disabled_state()

    def on_add_primer_row(self, idx: int) -> None:
        """Add a new empty primer row immediately below the row at idx.

        Syncs current UI to state, inserts a new empty primer dict,
        rebuilds the UI, and triggers the change handler.

        Args:
            idx: Zero-based index of the row below which to insert.
        """
        self.owner.sync_to_state(rebuild_if_needed=False)
        self.owner.input_data.primers.insert(
            idx + 1, {"name": "", "seq": "", "active": False}
        )
        self.owner.update_ui()
        if self.owner.on_change_handler is not None:
            self.owner.on_change_handler(None)

    def move_primers(self, indices: set[int], direction: int) -> None:
        """Move multiple primers at indices up or down.

        Args:
            indices: Set of zero-based indices of the primers to move.
            direction: -1 to move up, 1 to move down.
        """
        self.owner.sync_to_state(rebuild_if_needed=False)
        primers = self.owner.input_data.primers

        if direction == -1:
            sorted_sel = sorted(indices)
            if sorted_sel and sorted_sel[0] == 0:
                return
            for idx in sorted_sel:
                primers[idx], primers[idx - 1] = primers[idx - 1], primers[idx]
            self.owner.selected_indices = {i - 1 for i in indices}
            if self.owner.focused_primer_index is not None:
                self.owner.focused_primer_index -= 1
        elif direction == 1:
            sorted_sel = sorted(indices, reverse=True)
            if sorted_sel and sorted_sel[0] == len(primers) - 1:
                return
            for idx in sorted_sel:
                primers[idx], primers[idx + 1] = primers[idx + 1], primers[idx]
            self.owner.selected_indices = {i + 1 for i in indices}
            if self.owner.focused_primer_index is not None:
                self.owner.focused_primer_index += 1

        self.owner.update_ui()
        if self.owner.on_change_handler is not None:
            self.owner.on_change_handler(None)

    def delete_primers(self, indices_to_delete: set[int]) -> None:
        """Delete primers at indices and re-index controls in-place.

        Removes the specified primers from state and UI controls,
        adjusts the focused index, and re-indexes remaining rows.

        Args:
            indices_to_delete: Set of zero-based indices to remove.
        """
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
        min_deleted = min(indices_to_delete)
        new_len = len(new_primers)
        if new_len == 0 or (
            new_len == 1
            and not new_primers[0].get("name")
            and not new_primers[0].get("seq")
        ):
            self.owner.focused_primer_index = None
            self.owner.selected_indices = set()
        else:
            new_focus = min(min_deleted, new_len - 1)
            self.owner.focused_primer_index = new_focus
            self.owner.selected_indices = {new_focus}

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
                row.update_index(
                    new_idx=new_idx,
                    on_row_click=self.handle_row_click,
                )

        self.owner._update_row_highlights()
        self.owner._update_primer_info_panel()
        self.owner._update_header_checkbox_state()
        self.owner._update_delete_button_disabled_state()
        if self.owner.app_page:
            self.owner.app_page.update()

        if self.owner.on_change_handler is not None:
            self.owner.on_change_handler(None)

    def handle_drag_start(self, start_idx: int, e: ft.DragStartEvent) -> None:
        """Handle start of drag selection on a primer row."""
        self.drag_start_index = start_idx
        self.current_drag_y = 0.0
        ctrl = getattr(self.owner, "ctrl_pressed", False)
        if not ctrl:
            self.owner.selected_indices = {start_idx}
        else:
            self.owner.selected_indices.add(start_idx)
        self.owner.focused_primer_index = start_idx
        self.owner._update_row_highlights()
        self.owner._update_delete_button_disabled_state()

    def handle_drag_update(self, start_idx: int, e: ft.DragUpdateEvent) -> None:
        """Handle drag selection update over primer rows."""
        from .primer_row import PrimerRow

        intervals = []
        current_y = 0.0
        for row in self.owner.primers_list.controls:
            if isinstance(row, PrimerRow):
                h = (
                    55.0
                    if (row.name_field.error or row.seq_field.error)
                    else 30.0
                )
                intervals.append((current_y, current_y + h))
                current_y += h

        if not intervals or start_idx >= len(intervals):
            return

        delta_y = (
            getattr(e.local_delta, "y", 0.0)
            if getattr(e, "local_delta", None)
            else 0.0
        )
        self.current_drag_y += delta_y

        start_top_y = intervals[start_idx][0]
        start_bottom_y = intervals[start_idx][1]
        start_center_y = (start_top_y + start_bottom_y) / 2.0

        current_y_pos = start_center_y + self.current_drag_y

        target_idx = start_idx
        for idx, (top, bottom) in enumerate(intervals):
            if top <= current_y_pos < bottom:
                target_idx = idx
                break
        else:
            if current_y_pos < 0:
                target_idx = 0
            elif current_y_pos >= intervals[-1][1]:
                target_idx = len(intervals) - 1

        start = min(start_idx, target_idx)
        end = max(start_idx, target_idx)

        ctrl = getattr(self.owner, "ctrl_pressed", False)
        if not ctrl:
            self.owner.selected_indices = set(range(start, end + 1))
        else:
            for i in range(start, end + 1):
                self.owner.selected_indices.add(i)

        self.owner._update_row_highlights()
        self.owner._update_delete_button_disabled_state()

    def handle_drag_end(self, start_idx: int, e: ft.DragEndEvent) -> None:
        """Handle end of drag selection on a primer row."""
        self.owner._update_row_highlights()
        self.owner._update_primer_info_panel()
        self.owner._update_delete_button_disabled_state()
        if self.owner.app_page:
            self.owner.app_page.update()
