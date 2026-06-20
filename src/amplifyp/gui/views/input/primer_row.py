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

"""A single row representing a primer in the list."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, cast

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings

if TYPE_CHECKING:
    pass


class PrimerRow(ft.Container):  # type: ignore[misc]
    """A single row representing a primer in the list."""

    def __init__(
        self,
        idx: int,
        name: str,
        seq: str,
        is_active: bool,
        is_dup: bool,
        name_error: str | None,
        seq_error: str | None,
        font_family: str,
        name_column_width: float,
        settings: GUISettings,
        on_change_handler: Callable[[ft.Event | None], None],
        handle_field_focus: Callable[[ft.ControlEvent], None],
        handle_field_blur: Callable[[ft.ControlEvent], None],
        handle_field_submit: Callable[[ft.Event], None],
        on_row_click: Callable[[int, ft.TextField], None],
        on_move_primer: Callable[[int, int], None],
        on_delete_primer: Callable[[int], None],
        on_add_primer: Callable[[int], None],
        on_divider_pan: Callable[[ft.DragUpdateEvent], None],
        on_divider_pan_end: Callable[[ft.DragEndEvent], None],
        is_focused: bool,
        is_last_row: bool,
    ) -> None:
        """Initialise the PrimerRow.

        Args:
            idx: The zero-based index of this primer row.
            name: The primer name.
            seq: The primer DNA sequence.
            is_active: Whether this primer is selected/active.
            is_dup: Whether this primer is a duplicate.
            name_error: Error message for the name field, or None.
            seq_error: Error message for the sequence field, or None.
            font_family: Font family for sequence display.
            name_column_width: Width of the name column in pixels.
            settings: Application GUI settings instance.
            on_change_handler: Callback for field change events.
            handle_field_focus: Callback for field focus events.
            handle_field_blur: Callback for field blur events.
            handle_field_submit: Callback for field submit events.
            on_row_click: Callback when the row container is clicked.
            on_move_primer: Callback to move a primer up or down.
            on_delete_primer: Callback to delete a primer.
            on_add_primer: Callback to add a new primer below.
            on_divider_pan: Callback for dragging the name/sequence divider.
            on_divider_pan_end: Callback for ending the divider drag.
            is_focused: Whether this row is currently focused.
            is_last_row: Whether this is the last row in the list.
        """
        has_err = bool(name_error or seq_error)
        super().__init__(
            data=idx,
            bgcolor=GUIColours.DUPLICATE_BG if is_dup else None,
            padding=0,
            height=30 if not has_err else None,
        )
        self.idx = idx
        self.settings = settings
        show_temp = self.settings.get("show_primer_temperature", False)

        tm_val = ""
        self._tm_value: float | None = None
        if seq.strip():
            try:
                from amplifyp.dna import Primer
                from amplifyp.gui.util import clean_sequence

                cleaned_seq = clean_sequence(seq)
                if cleaned_seq:
                    primer_obj = Primer(sequence=cleaned_seq, name=name)
                    tm = self.settings.calculate_primer_tm(primer_obj)
                    self._tm_value = tm
                    tm_val = f"{tm:.1f}°C"
            except Exception:
                tm_val = "-"

        scheme = self.settings.get("tm_colour_scheme", "None")
        from amplifyp.gui.colours import tm_colour

        _tm_colour = (
            tm_colour(self._tm_value, scheme)
            if self._tm_value is not None
            else None
        )
        self.tm_text = ft.Text(
            value=tm_val,
            size=self.settings.get("font_size_body", 13),
            color=_tm_colour,
            selectable=False,
        )
        self.tm_container = ft.Container(
            content=self.tm_text,
            width=50,
            padding=ft.Padding(5, 0, 0, 0),
            alignment=ft.Alignment(-1, 0),
            visible=show_temp,
        )
        self.tm_divider = ft.Container(
            width=4,
            bgcolor=GUIColours.DIVIDER_GREY,
            margin=0,
            height=30,
            visible=show_temp,
        )

        is_empty = not name.strip() or not seq.strip()
        self.checkbox = ft.Checkbox(
            value=is_active if not is_empty else False,
            on_change=on_change_handler,
            disabled=is_empty,
            visible=True,
        )
        self.checkbox_container = ft.Container(
            content=self.checkbox,
            width=55,
            height=30,
            alignment=ft.Alignment(0, 0),
        )
        self.name_field = ft.TextField(
            value=name,
            hint_text="New Primer Name",
            dense=True,
            content_padding=ft.Padding(5, 0, 0, 0),
            height=30 if not name_error else None,
            width=name_column_width,
            border=ft.InputBorder.NONE,
            data={"idx": idx, "field": "name"},
            on_focus=handle_field_focus,
            on_blur=handle_field_blur,
            on_submit=handle_field_submit,
        )
        self.seq_field = ft.TextField(
            value=seq,
            hint_text="New Primer Sequence",
            dense=True,
            content_padding=ft.Padding(5, 0, 5, 0),
            height=30 if not seq_error else None,
            border=ft.InputBorder.NONE,
            text_style=ft.TextStyle(font_family=font_family),
            data={"idx": idx, "field": "seq"},
            expand=True,
            on_focus=handle_field_focus,
            on_blur=handle_field_blur,
            on_submit=handle_field_submit,
        )
        if name_error:
            self.name_field.error = name_error
        if seq_error:
            self.seq_field.error = seq_error

        self.divider = ft.GestureDetector(
            on_pan_update=on_divider_pan,
            on_pan_end=on_divider_pan_end,
            content=ft.Container(
                width=4,
                bgcolor=GUIColours.DIVIDER_GREY,
                margin=0,
                height=30,
            ),
            mouse_cursor=ft.MouseCursor.RESIZE_LEFT_RIGHT,
        )

        self.active_divider = ft.Container(
            width=4,
            bgcolor=GUIColours.DIVIDER_GREY,
            margin=0,
            height=30,
        )

        self.reorder_controls = None
        self.control_container = None
        up_button = ft.IconButton(
            icon=ft.Icons.ARROW_UPWARD,
            icon_size=16,
            width=24,
            height=24,
            padding=0,
            tooltip="Move Up",
            disabled=(idx == 0),
            on_click=lambda e: on_move_primer(idx, -1),
        )
        down_button = ft.IconButton(
            icon=ft.Icons.ARROW_DOWNWARD,
            icon_size=16,
            width=24,
            height=24,
            padding=0,
            tooltip="Move Down",
            disabled=is_last_row,
            on_click=lambda e: on_move_primer(idx, 1),
        )
        delete_button = ft.IconButton(
            icon=ft.Icons.DELETE_OUTLINE,
            icon_size=16,
            width=24,
            height=24,
            padding=0,
            tooltip="Delete Primer",
            on_click=lambda e: on_delete_primer(idx),
        )
        add_button = ft.IconButton(
            icon=ft.Icons.ADD_CIRCLE_OUTLINE,
            icon_size=16,
            width=24,
            height=24,
            padding=0,
            tooltip="Add Primer Below",
            on_click=lambda e: on_add_primer(idx),
        )
        self.reorder_controls = ft.Row(
            [add_button, delete_button, up_button, down_button],
            spacing=2,
            alignment=ft.MainAxisAlignment.CENTER,
        )
        self.control_container = ft.Container(
            content=self.reorder_controls,
            width=108 if is_focused else 0,
            height=30,
            alignment=ft.Alignment(0, 0),
        )
        self.reorder_controls.visible = is_focused

        controls = [
            self.checkbox_container,
            self.active_divider,
            self.name_field,
            self.divider,
        ]
        if show_temp:
            controls.extend([self.tm_container, self.tm_divider])
        controls.extend([self.seq_field, self.control_container])

        self.content = ft.Row(
            controls,
            spacing=0,
            vertical_alignment=ft.CrossAxisAlignment.START,
        )
        self.on_click = lambda e: on_row_click(idx, self.name_field)

    def update_highlight_and_reorder(
        self, is_focused: bool, is_dup: bool
    ) -> None:
        """Update the background colour and reorder buttons layout.

        Args:
            is_focused: Whether this row should show the reorder controls.
            is_dup: Whether this primer is a duplicate.
        """
        if is_focused:
            self.bgcolor = GUIColours.SELECTED_ROW_BG
        elif is_dup:
            self.bgcolor = GUIColours.DUPLICATE_BG
        else:
            self.bgcolor = None  # type: ignore[assignment]

        if self.control_container is not None:
            self.control_container.width = 108 if is_focused else 0
            try:
                self.control_container.update()
            except RuntimeError:
                pass
        if self.reorder_controls is not None:
            self.reorder_controls.visible = is_focused
            try:
                self.reorder_controls.update()
            except RuntimeError:
                pass

    def set_error(self, err: dict[str, str | None] | str | None) -> None:
        """Set or clear the error message.

        Also adjusts the sequence field, name field and container height.

        Args:
            err: Error dict with 'name' and 'seq' keys, a string error
                message, or None to clear errors.
        """
        if isinstance(err, dict):
            name_error = err.get("name")
            seq_error = err.get("seq")
        else:
            # Compatibility: if string or None is passed
            if err == "Duplicate primer name":
                name_error = err
                seq_error = None
            else:
                name_error = None
                seq_error = err

        self.name_field.error = name_error
        self.name_field.height = 30 if not name_error else None
        self.seq_field.error = seq_error
        self.seq_field.height = 30 if not seq_error else None

        has_err = bool(name_error or seq_error)
        self.height = 30 if not has_err else None

        is_empty = (
            not self.name_field.value.strip()
            or not self.seq_field.value.strip()
        )
        self.checkbox.disabled = is_empty
        if self.checkbox.disabled:
            self.checkbox.value = False

    def update_index(
        self,
        new_idx: int,
        is_last_row: bool,
        on_move_primer: Callable[[int, int], None],
        on_delete_primer: Callable[[int], None],
        on_add_primer: Callable[[int], None],
        on_row_click: Callable[[int, ft.TextField], None],
    ) -> None:
        """Update the index of the row and refresh its handlers and controls.

        Args:
            new_idx: The new zero-based index for this primer row.
            is_last_row: Whether this row is now the last in the list.
            on_move_primer: Callback to move a primer up or down.
            on_delete_primer: Callback to delete a primer.
            on_add_primer: Callback to add a new primer below.
            on_row_click: Callback when the row container is clicked.
        """
        self.data = new_idx
        self.idx = new_idx
        self.name_field.data = {"idx": new_idx, "field": "name"}
        self.seq_field.data = {"idx": new_idx, "field": "seq"}

        # Update click handler with the new index
        self.on_click = lambda e: on_row_click(new_idx, self.name_field)

        # Update reorder control buttons with the new index and state
        if self.reorder_controls is not None:
            add_button = cast(ft.IconButton, self.reorder_controls.controls[0])
            delete_button = cast(
                ft.IconButton, self.reorder_controls.controls[1]
            )
            up_button = cast(ft.IconButton, self.reorder_controls.controls[2])
            down_button = cast(ft.IconButton, self.reorder_controls.controls[3])

            add_button.on_click = lambda e: on_add_primer(new_idx)
            delete_button.on_click = lambda e: on_delete_primer(new_idx)
            up_button.on_click = lambda e: on_move_primer(new_idx, -1)
            up_button.disabled = new_idx == 0
            down_button.on_click = lambda e: on_move_primer(new_idx, 1)
            down_button.disabled = is_last_row

    def update_tm(self, settings: GUISettings) -> None:
        """Update the displayed Tm in-place based on the current sequence."""
        # ponytail: calculates Tm using user-selected formula
        seq_val = self.seq_field.value
        name_val = self.name_field.value
        tm_val = ""
        self._tm_value = None
        if seq_val and seq_val.strip():
            try:
                from amplifyp.dna import Primer
                from amplifyp.gui.util import clean_sequence

                cleaned_seq = clean_sequence(seq_val)
                if cleaned_seq:
                    primer_obj = Primer(sequence=cleaned_seq, name=name_val)
                    tm = settings.calculate_primer_tm(primer_obj)
                    self._tm_value = tm
                    tm_val = f"{tm:.1f}°C"
            except Exception:
                tm_val = "-"
        self.tm_text.value = tm_val
        scheme = settings.get("tm_colour_scheme", "None")
        from amplifyp.gui.colours import tm_colour

        self.tm_text.color = (
            tm_colour(self._tm_value, scheme)
            if self._tm_value is not None
            else None
        )
        try:
            self.tm_text.update()
        except RuntimeError:
            pass
