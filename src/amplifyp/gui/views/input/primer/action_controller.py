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

from __future__ import annotations

import inspect
import logging
from typing import TYPE_CHECKING

import flet as ft

from amplifyp.gui.utils.gui_helpers import focus_async

from .row import PrimerRow

if TYPE_CHECKING:
    from .input import PrimerInput

logger = logging.getLogger(__name__)


class PrimerActionController:
    """Handles mutation operations on the primers list."""

    def __init__(self, owner: PrimerInput) -> None:
        """Initialise the PrimerActionController.

        Args:
            owner: The parent PrimerInput component that owns this controller.
        """
        self.owner = owner
        self.current_drag_y = 0.0
        self._click_a: int | None = None
        self._click_b: int | None = None

    def handle_row_click(self, idx: int, name_edit: ft.TextField) -> None:
        """Handle single click on the row container.

        Immediately toggles selection of the row and resets
        click_a/click_b, clearing any pending double-click range.

        Args:
            idx: Zero-based index of the clicked primer row.
            name_edit: The name TextField control to focus.
        """
        owner = self.owner
        self._click_a = None
        self._click_b = None

        if idx in owner.selected_indices:
            owner.selected_indices.discard(idx)
            if owner.focused_primer_index == idx:
                owner.focused_primer_index = None
        else:
            owner.selected_indices.add(idx)
            owner.focused_primer_index = idx

            res = name_edit.focus()
            if inspect.iscoroutine(res) and owner.app_page:
                owner.app_page.run_task(focus_async, res)

        owner._update_row_highlights()
        owner._update_primer_info_panel()
        owner._update_delete_button_disabled_state()

    def handle_row_double_click(
        self, idx: int, _name_edit: ft.TextField
    ) -> None:
        """Handle double click on the row container.

        On the first double-click, toggles selection of the clicked row
        and sets click_a. On the second double-click, sets click_b,
        toggles highlighting between click_a (exclusive) and click_b
        (inclusive), then resets both to None.

        Args:
            idx: Zero-based index of the clicked primer row.
            _name_edit: The name TextField control (unused).
        """
        owner = self.owner

        # Reset anchor if it is out of bounds due to list mutations
        if self._click_a is not None and not (
            0 <= self._click_a < len(owner.input_data.primers)
        ):
            self._click_a = None

        if self._click_a is None:
            # First double-click: toggle the row, set click_a
            if idx in owner.selected_indices:
                owner.selected_indices.discard(idx)
                if owner.focused_primer_index == idx:
                    owner.focused_primer_index = None
            else:
                owner.selected_indices.add(idx)
                owner.focused_primer_index = idx
            self._click_a = idx
        else:
            # Second double-click: set click_b, toggle range
            # exclusive of click_a (already selected)
            self._click_b = idx
            click_a = self._click_a
            start = min(click_a, self._click_b)
            end = max(click_a, self._click_b)
            for i in range(start, end + 1):
                if i == click_a:
                    continue
                if i in owner.selected_indices:
                    owner.selected_indices.discard(i)
                else:
                    owner.selected_indices.add(i)
            self._click_a = None
            self._click_b = None

        owner._update_row_highlights()
        owner._update_primer_info_panel()
        owner._update_delete_button_disabled_state()

    def on_add_primer_row(self, idx: int) -> None:
        """Add a new empty primer row immediately below the row at idx.

        Syncs current UI to state, inserts a new empty primer dict,
        rebuilds the UI, and triggers the change handler.

        Args:
            idx: Zero-based index of the row below which to insert.
        """
        self._click_a = None
        self._click_b = None
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
        self._click_a = None
        self._click_b = None
        self.owner.sync_to_state(rebuild_if_needed=False)
        primers = self.owner.input_data.primers
        valid_indices = {i for i in indices if 0 <= i < len(primers)}
        if not valid_indices:
            return

        if direction == -1:
            sorted_sel = sorted(valid_indices)
            if sorted_sel and sorted_sel[0] == 0:
                return
            for idx in sorted_sel:
                primers[idx], primers[idx - 1] = primers[idx - 1], primers[idx]
            self.owner.selected_indices = {i - 1 for i in valid_indices}
            if self.owner.focused_primer_index in valid_indices:
                self.owner.focused_primer_index -= 1
        elif direction == 1:
            sorted_sel = sorted(valid_indices, reverse=True)
            if sorted_sel and sorted_sel[0] == len(primers) - 1:
                return
            for idx in sorted_sel:
                primers[idx], primers[idx + 1] = primers[idx + 1], primers[idx]
            self.owner.selected_indices = {i + 1 for i in valid_indices}
            if self.owner.focused_primer_index in valid_indices:
                self.owner.focused_primer_index += 1

        self.owner.update_ui()
        if self.owner.on_change_handler is not None:
            self.owner.on_change_handler(None)

    def delete_primers(self, indices_to_delete: set[int]) -> None:
        """Delete primers asynchronously to avoid Flet race conditions."""
        self.owner.sync_to_state(rebuild_if_needed=False)
        primers = self.owner.input_data.primers
        primers_to_delete = {
            id(primers[i]) for i in indices_to_delete if 0 <= i < len(primers)
        }

        page = None
        try:
            page = self.owner.page
        except RuntimeError as e:
            logger.debug("Failed to get owner page: %s", e)
        if page:

            async def delayed_delete() -> None:
                import asyncio

                await asyncio.sleep(0.05)
                self._delete_primers_impl(primers_to_delete)

            page.run_task(delayed_delete)
        else:
            self._delete_primers_impl(primers_to_delete)

    def _delete_primers_impl(self, primers_to_delete: set[int]) -> None:
        """Perform the actual deletion of primers from the input data.

        This method removes the primers identified by their object IDs, updates
        the focused primer index, and re-indexes the remaining primer rows
        in the UI.
        """
        self._click_a = None
        self._click_b = None
        if not primers_to_delete:
            return

        primers = self.owner.input_data.primers
        deleted_indices = {
            i for i, p in enumerate(primers) if id(p) in primers_to_delete
        }
        if not deleted_indices:
            return

        # Keep only primers NOT in the deleted set
        new_primers = [p for p in primers if id(p) not in primers_to_delete]
        if not new_primers:
            new_primers = [{"name": "", "seq": "", "active": False}]
        self.owner.input_data.primers = new_primers

        # Adjust focus index
        min_deleted = min(deleted_indices)
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

        # Filter the controls list
        remaining_controls = []
        for row in self.owner.primers_list.controls:
            if isinstance(row, PrimerRow):
                if row.data in deleted_indices:
                    continue
                remaining_controls.append(row)

        if not remaining_controls:
            self.owner.primers_list.controls.clear()
            self.owner.update_ui()
            return

        self.owner.primers_list.controls = remaining_controls
        num_remaining = len(remaining_controls)
        # Find the earliest deleted index to start re-indexing from
        start_reindex = min(deleted_indices)

        for new_idx in range(start_reindex, num_remaining):
            row = remaining_controls[new_idx]
            if isinstance(row, PrimerRow):
                row.update_index(
                    new_idx=new_idx,
                    on_row_click=self.handle_row_click,
                    on_row_double_click=self.handle_row_double_click,
                )

        self.owner._update_row_highlights()
        self.owner._update_primer_info_panel()
        self.owner._update_header_checkbox_state()
        self.owner._update_delete_button_disabled_state()
        if self.owner.on_change_handler is not None:
            self.owner.on_change_handler(None)
        elif self.owner.app_page:
            self.owner.app_page.update()

    def handle_drag_start(self, start_idx: int, _e: ft.DragStartEvent) -> None:
        """Handle start of drag reordering on a primer row."""
        primers = self.owner.input_data.primers
        self.owner.selected_indices = {
            i for i in self.owner.selected_indices if 0 <= i < len(primers)
        }
        if start_idx not in self.owner.selected_indices:
            self.owner.selected_indices = {start_idx}
            self.owner.focused_primer_index = start_idx

        # Find the contiguous block of highlighted rows containing start_idx
        selected = sorted(self.owner.selected_indices)
        block_start = start_idx
        block_end = start_idx

        idx_pos = selected.index(start_idx)
        for i in range(idx_pos - 1, -1, -1):
            if selected[i] == block_start - 1:
                block_start = selected[i]
            else:
                break

        for i in range(idx_pos + 1, len(selected)):
            if selected[i] == block_end + 1:
                block_end = selected[i]
            else:
                break

        self.drag_block = list(range(block_start, block_end + 1))
        self.current_drag_y = 0.0

        self.owner._update_row_highlights()
        self.owner._update_delete_button_disabled_state()

    def handle_drag_update(
        self, _start_idx: int, e: ft.DragUpdateEvent
    ) -> None:
        """Handle live drag-and-drop reordering of contiguous rows."""
        self._click_a = None
        self._click_b = None
        if not hasattr(self, "drag_block") or not self.drag_block:
            return

        delta_y = (
            getattr(e.local_delta, "y", 0.0)
            if getattr(e, "local_delta", None)
            else getattr(e, "delta_y", 0.0)
        )
        self.current_drag_y += delta_y

        primers = self.owner.input_data.primers
        changed = False

        controls = self.owner.primers_list.controls

        def get_row_height(idx: int) -> float:
            """Determine the height of a primer row.

            Args:
                idx: The index of the row.

            Returns:
                The height in pixels (55.0 for error rows, 30.0 otherwise).
            """
            if 0 <= idx < len(controls):
                row = controls[idx]
                if isinstance(row, PrimerRow):
                    return (
                        55.0
                        if (row.name_field.error or row.seq_field.error)
                        else 30.0
                    )
            return 30.0

        # Try to move block down
        while True:
            block_end = self.drag_block[-1]
            if block_end >= len(primers) - 1:
                break

            row_below_height = get_row_height(block_end + 1)
            if self.current_drag_y > row_below_height / 2.0:
                target_idx = block_end + 1
                if target_idx >= len(primers) or self.drag_block[0] >= len(
                    primers
                ):
                    break
                # Move the row below the block to above the block
                primers.insert(self.drag_block[0], primers.pop(target_idx))

                self.drag_block = [i + 1 for i in self.drag_block]
                self.current_drag_y -= row_below_height
                changed = True
            else:
                break

        # Try to move block up
        while True:
            block_start = self.drag_block[0]
            if block_start <= 0 or block_start >= len(primers):
                break

            row_above_height = get_row_height(block_start - 1)
            if self.current_drag_y < -row_above_height / 2.0:
                target_idx = block_start - 1
                if target_idx >= len(primers) or self.drag_block[-1] >= len(
                    primers
                ):
                    break
                # Move the row above the block to below the block
                primers.insert(self.drag_block[-1], primers.pop(target_idx))

                self.drag_block = [i - 1 for i in self.drag_block]
                self.current_drag_y += row_above_height
                changed = True
            else:
                break

        if changed:
            self.owner.selected_indices = set(self.drag_block)
            if self.owner.focused_primer_index is not None:
                self.owner.focused_primer_index = self.drag_block[0]
            self.owner.update_ui()

    def handle_drag_end(self, _start_idx: int, _e: ft.DragEndEvent) -> None:
        """Handle end of drag reordering on a primer row."""
        if hasattr(self, "drag_block"):
            delattr(self, "drag_block")

        self.owner._update_row_highlights()
        self.owner._update_primer_info_panel()
        self.owner._update_delete_button_disabled_state()
        if self.owner.on_change_handler is not None:
            self.owner.on_change_handler(None)
        elif self.owner.app_page:
            self.owner.app_page.update()

    def header_add_click(self, e: ft.Event | None) -> None:
        """Handle header Add button click."""
        num_primers = len(self.owner.input_data.primers)
        if self.owner.selected_indices:
            idx = min(max(self.owner.selected_indices), num_primers - 1)
        elif self.owner.focused_primer_index is not None:
            idx = min(self.owner.focused_primer_index, num_primers - 1)
        else:
            idx = num_primers - 1

        self.owner.selected_indices = {idx + 1}
        self.owner.focused_primer_index = idx + 1
        self.on_add_primer_row(idx)

    def header_delete_click(self, e: ft.Event | None) -> None:
        """Handle header Delete button click."""
        if self.owner.selected_indices:
            self.delete_primers(self.owner.selected_indices.copy())
            self.owner._update_header_buttons_state()

    def header_up_click(self, e: ft.Event | None) -> None:
        """Handle header Move Up button click."""
        if self.owner.selected_indices and min(self.owner.selected_indices) > 0:
            self.move_primers(self.owner.selected_indices, -1)
            self.owner._update_header_buttons_state()

    def header_down_click(self, e: ft.Event | None) -> None:
        """Handle header Move Down button click."""
        if (
            self.owner.selected_indices
            and max(self.owner.selected_indices)
            < len(self.owner.input_data.primers) - 1
        ):
            self.move_primers(self.owner.selected_indices, 1)
            self.owner._update_header_buttons_state()
