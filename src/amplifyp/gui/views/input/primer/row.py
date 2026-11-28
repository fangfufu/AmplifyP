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

"""A single row representing a primer in the list."""

from __future__ import annotations

import asyncio
import inspect
import logging
from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.dna import Primer
from amplifyp.gui.colours import GUIColours, tm_colour
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.data_helpers import clean_sequence

logger = logging.getLogger(__name__)


def _format_variant_count(variant_count: int) -> str:
    """Format an arbitrary-size variant count without float overflow."""
    if variant_count < 1_000_000:
        return f"{variant_count:,}"
    digits = str(variant_count)
    exp = len(digits) - 1
    mantissa = round(int(digits[:4]) / 1000.0, 2)
    if mantissa >= 10.0:
        mantissa /= 10.0
        exp += 1
    return f"{mantissa:.2f}e+{exp:02d}"


class PrimerRow(ft.Container):  # type: ignore[misc]
    """A single row representing a primer in the list."""

    idx: int
    name_field: ft.TextField
    seq_field: ft.TextField

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
        handle_field_focus: Callable[[ft.Event[ft.TextField]], None],
        handle_field_blur: Callable[[ft.Event[ft.TextField]], None],
        handle_field_submit: Callable[[ft.Event[ft.TextField]], None],
        on_row_click: Callable[[int, ft.TextField], None],
        on_row_double_click: Callable[[int, ft.TextField], None],
        on_divider_pan: Callable[[ft.DragUpdateEvent], None],
        on_divider_pan_end: Callable[[ft.DragEndEvent], None],
        is_focused: bool,
        is_last_row: bool,
        on_drag_start: Callable[[int, ft.DragStartEvent], None] | None = None,
        on_drag_update: Callable[[int, ft.DragUpdateEvent], None] | None = None,
        on_drag_end: Callable[[int, ft.DragEndEvent], None] | None = None,
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
            on_row_click: Callback when the row container is single-clicked.
            on_row_double_click: Callback for double-click on row.
            on_divider_pan: Callback for dragging the name/sequence divider.
            on_divider_pan_end: Callback for ending the divider drag.
            is_focused: Whether this row is currently focused.
            is_last_row: Whether this is the last row in the list.
            on_drag_start: Callback for starting a drag-selection.
            on_drag_update: Callback for updating a drag-selection.
            on_drag_end: Callback for ending a drag-selection.
        """
        init_msgs = [m for m in (name_error, seq_error) if m]
        init_tooltip = "; ".join(init_msgs) if init_msgs else None
        super().__init__(
            data=idx,
            bgcolor=GUIColours.DUPLICATE_BG if is_dup else None,
            padding=0,
            height=30,
            tooltip=init_tooltip,
        )
        self.idx = idx
        self.settings = settings
        self._handle_field_blur = handle_field_blur
        self.is_last_row = is_last_row
        show_temp = self.settings.get("show_primer_temperature", False)
        scheme = self.settings.get("tm_colour_scheme", "None")
        self._last_show_temp = show_temp
        self._tm_value, tm_val, tm_tooltip = self._compute_tm_display(
            seq, name, show_temp, self.settings
        )
        self._last_tm_cache_key = self._tm_cache_key(
            seq, name, show_temp, scheme, self.settings
        )

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
            padding=ft.Padding(0, 0, 5, 0),
            alignment=ft.Alignment(1, 0),
            visible=show_temp,
            tooltip=tm_tooltip,
        )
        self.tm_divider = ft.Container(
            width=4,
            bgcolor=GUIColours.DIVIDER_GREY,
            margin=0,
            height=30,
            visible=show_temp,
        )

        self.checkbox = ft.Checkbox(
            value=is_active,
            on_change=on_change_handler,
        )
        self.checkbox_container = ft.Container(
            content=self.checkbox,
            alignment=ft.Alignment(0, 0),
            width=30,
            padding=0,
        )

        self._drag_icon = ft.Icon(
            ft.Icons.DRAG_INDICATOR,
            size=16,
            color=GUIColours.DIVIDER_GREY,
        )
        self.drag_handle = ft.GestureDetector(
            on_pan_start=lambda e: (
                on_drag_start(self.idx, e) if on_drag_start else None
            ),
            on_pan_update=lambda e: (
                on_drag_update(self.idx, e) if on_drag_update else None
            ),
            on_pan_end=lambda e: (
                on_drag_end(self.idx, e) if on_drag_end else None
            ),
            on_tap=lambda e: on_row_click(self.idx, self.name_field),
            on_double_tap=lambda e: on_row_double_click(
                self.idx,
                self.name_field,
            ),
            mouse_cursor=ft.MouseCursor.CLICK,
            content=ft.Container(
                content=self._drag_icon,
                alignment=ft.Alignment(0, 0),
                width=25,
                padding=0,
                height=30,
            ),
        )

        self.name_field = ft.TextField(
            value=name,
            dense=True,
            content_padding=ft.Padding(5, 0, 0, 0),
            height=30,
            border=ft.NoInputBorder(),
            multiline=True,
            ignore_up_down_keys=True,
            fit_parent_size=True,
            data={"idx": idx, "field": "name", "cursor_pos": 0},
            on_focus=handle_field_focus,
            on_blur=handle_field_blur,
            on_submit=handle_field_submit,
            on_change=on_change_handler,
        )
        self.name_container = ft.Container(
            content=self.name_field,
            width=1000,
            height=30,
        )
        self.name_scroll = ft.ListView(
            horizontal=True,
            width=name_column_width,
            height=30,
            controls=[self.name_container],
        )

        self.seq_field = ft.TextField(
            value=seq,
            dense=True,
            content_padding=ft.Padding(5, 0, 5, 0),
            height=30,
            border=ft.NoInputBorder(),
            text_style=ft.TextStyle(font_family=font_family),
            multiline=True,
            ignore_up_down_keys=True,
            fit_parent_size=True,
            data={"idx": idx, "field": "seq", "cursor_pos": 0},
            on_focus=handle_field_focus,
            on_blur=handle_field_blur,
            on_submit=handle_field_submit,
            on_change=on_change_handler,
            selection=ft.TextSelection(base_offset=0, extent_offset=0),
        )
        self.seq_container = ft.Container(
            content=self.seq_field,
            width=5000,
            height=30,
        )
        self.seq_scroll = ft.ListView(
            horizontal=True,
            expand=True,
            height=30,
            controls=[self.seq_container],
        )

        self.name_field.on_selection_change = self._on_selection_change
        self.seq_field.on_selection_change = self._on_selection_change

        self.name_field.on_blur = self._on_blur
        self.seq_field.on_blur = self._on_blur

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

        row_body = [
            self.checkbox_container,
            self.active_divider,
            self.name_scroll,
            self.divider,
            self.seq_scroll,
        ]
        if show_temp:
            row_body.extend([self.tm_divider, self.tm_container])

        self._row_content = ft.Row(
            row_body,
            spacing=0,
            vertical_alignment=ft.CrossAxisAlignment.START,
        )

        self._row_gesture_detector = ft.GestureDetector(
            on_tap=lambda e: on_row_click(self.idx, self.name_field),
            mouse_cursor=ft.MouseCursor.CLICK,
            content=self._row_content,
        )

        self.content = ft.Row(
            [self.drag_handle, self._row_gesture_detector],
            spacing=0,
            vertical_alignment=ft.CrossAxisAlignment.START,
        )
        self.height = 30

    def _on_selection_change(self, e: ft.Event) -> None:
        """Handle selection change for name or seq field."""
        sel = (
            getattr(e.control, "selection", None)
            if hasattr(e, "control")
            else None
        )
        if sel is None or not hasattr(sel, "base_offset"):
            return
        field_data: dict[str, str] = e.control.data
        field_data["cursor_pos"] = sel.base_offset

    def _on_blur(self, e: ft.Event[ft.TextField]) -> None:
        """Handle blur for name or seq field, scrolling to top first."""
        scroll_target: ft.ListView = (
            self.name_scroll
            if e.control.data["field"] == "name"
            else self.seq_scroll
        )
        res = scroll_target.scroll_to(offset=0)
        page = e.page or self.page
        if inspect.iscoroutine(res):

            async def _do_scroll() -> None:
                try:
                    await res
                except RuntimeError:
                    pass
                try:
                    if scroll_target.page:
                        scroll_target.update()
                except (RuntimeError, AssertionError):
                    pass

            scheduled = False
            if page and callable(getattr(page, "run_task", None)):
                try:
                    page.run_task(_do_scroll)  # pyright: ignore[reportAttributeAccessIssue]
                    scheduled = True
                except (RuntimeError, AttributeError):
                    pass
            if not scheduled:
                try:
                    _task = asyncio.get_running_loop().create_task(_do_scroll())
                    _task.add_done_callback(lambda t: t.exception())
                    scheduled = True
                except RuntimeError:
                    pass
            if not scheduled:
                _loop = asyncio.new_event_loop()
                try:
                    _loop.run_until_complete(_do_scroll())
                except RuntimeError:
                    pass
                finally:
                    _loop.close()
        try:
            if scroll_target.page:
                scroll_target.update()
        except (RuntimeError, AssertionError):
            pass
        self._handle_field_blur(e)

    def update_highlight_and_reorder(
        self, is_focused: bool, is_dup: bool, is_error: bool = False
    ) -> None:
        """Update the background colour.

        Args:
            is_focused: Whether this row is currently focused.
            is_dup: Whether this primer is a duplicate.
            is_error: Whether this primer has a validation error.
        """
        has_issue = is_dup or is_error
        if has_issue and is_focused:
            self.bgcolor = GUIColours.FOCUSED_ERROR_BG
        elif has_issue:
            self.bgcolor = GUIColours.ERROR_BG
        elif is_focused:
            self.bgcolor = GUIColours.SELECTED_ROW_BG
        else:
            self.bgcolor = None

    def set_error(self, err: dict[str, str | None] | str | None) -> None:
        """Update error state for the row.

        Field-level error messages are omitted in favour of row background
        colour highlights, tooltips, and error dialogues. Row height remains
        at 30px.

        Args:
            err: Error dict with 'name' and 'seq' keys, a string error
                message, or None to clear errors.
        """
        self.name_field.error = None
        self.seq_field.error = None
        if isinstance(err, dict):
            msgs = [m for m in (err.get("name"), err.get("seq")) if m]
            self.tooltip = "; ".join(msgs) if msgs else None
        else:
            self.tooltip = err or None
        self.checkbox.disabled = False

    def update_index(
        self,
        new_idx: int,
        on_row_click: Callable[[int, ft.TextField], None],
        on_row_double_click: Callable[[int, ft.TextField], None],
    ) -> None:
        """Update the index of the row and refresh its handlers and controls.

        Args:
            new_idx: The zero-based index for this primer row.
            on_row_click: Callback when the row container is clicked.
            on_row_double_click: Callback for double-click on row.
        """
        self.data = new_idx
        self.idx = new_idx
        name_cursor = self.name_field.data.get("cursor_pos", 0)
        seq_cursor = self.seq_field.data.get("cursor_pos", 0)
        self.name_field.data = {
            "idx": new_idx,
            "field": "name",
            "cursor_pos": name_cursor,
        }
        self.seq_field.data = {
            "idx": new_idx,
            "field": "seq",
            "cursor_pos": seq_cursor,
        }

        # Update click handler with the new index
        self._row_gesture_detector.on_tap = lambda e: on_row_click(
            new_idx, self.name_field
        )
        self._row_gesture_detector.on_double_tap = lambda e: (
            on_row_double_click(new_idx, self.name_field)
        )
        self.drag_handle.on_tap = lambda e: on_row_click(
            new_idx, self.name_field
        )
        self.drag_handle.on_double_tap = lambda e: on_row_double_click(
            new_idx, self.name_field
        )

    @staticmethod
    def _compute_tm_display(
        seq_val: str | None,
        name_val: str | None,
        show_temp: bool,
        settings: GUISettings,
    ) -> tuple[float | None, str, str | None]:
        """Compute raw, formatted, and tooltip Tm display values."""
        tm_value: float | None = None
        tm_val = ""
        tm_tooltip: str | None = None
        if show_temp and seq_val and seq_val.strip():
            try:
                cleaned_seq = clean_sequence(seq_val)
                if cleaned_seq:
                    primer_obj = Primer(
                        sequence=cleaned_seq, name=name_val or ""
                    )
                    tm_value, tm_range = settings.resolve_primer_tm_display(
                        primer_obj
                    )
                    tm_val = f"{tm_value:.1f}°C"
                    if tm_range is not None:
                        label = (
                            "approx. range"
                            if tm_range.is_approximate
                            else "range"
                        )
                        count_str = _format_variant_count(
                            tm_range.variant_count
                        )
                        tm_tooltip = (
                            f"Max Tm: {tm_value:.1f}°C ({label}: "
                            f"{tm_range.min_tm:.1f} to "
                            f"{tm_range.max_tm:.1f}°C, "
                            f"{count_str} variants)"
                        )
            except (
                ValueError,
                AttributeError,
                ArithmeticError,
                KeyError,
            ):
                logger.debug(
                    "Failed to calculate Tm for primer '%s'",
                    name_val,
                    exc_info=True,
                )
                tm_val = "-"
        return tm_value, tm_val, tm_tooltip

    @staticmethod
    def _tm_cache_key(
        seq_val: str | None,
        name_val: str | None,
        show_temp: bool,
        scheme: str,
        settings: GUISettings,
    ) -> tuple[Any, ...]:
        """Build a cache key covering primer inputs and Tm configuration."""
        if not show_temp:
            return (seq_val, name_val, False, scheme)
        return (
            seq_val,
            name_val,
            True,
            scheme,
            settings.get(
                "tm_method", "SantaLucia 1998 / Owczarzy 2008 (Default)"
            ),
            settings.get_tm_settings(),
        )

    def update_appearance(
        self, font_family: str, settings: GUISettings
    ) -> None:
        """Refresh sequence font family and theme divider colours."""
        self.settings = settings
        self.seq_field.text_style = ft.TextStyle(font_family=font_family)
        self._drag_icon.color = GUIColours.DIVIDER_GREY
        self.active_divider.bgcolor = GUIColours.DIVIDER_GREY
        if isinstance(self.divider.content, ft.Container):
            self.divider.content.bgcolor = GUIColours.DIVIDER_GREY
        self.tm_divider.bgcolor = GUIColours.DIVIDER_GREY
        self.tm_text.size = settings.get("font_size_body", 13)

    def update_tm(self, settings: GUISettings) -> None:
        """Update the displayed Tm in-place based on the current sequence."""
        self.settings = settings
        seq_val = self.seq_field.value
        name_val = self.name_field.value
        show_temp = settings.get("show_primer_temperature", False)
        scheme = settings.get("tm_colour_scheme", "None")

        # Check cache before recomputing
        cache_key = self._tm_cache_key(
            seq_val, name_val, show_temp, scheme, settings
        )
        if getattr(self, "_last_tm_cache_key", None) == cache_key:
            tm_val_stored = getattr(self, "_tm_value", None)
            self.tm_text.color = (
                tm_colour(tm_val_stored, scheme)
                if tm_val_stored is not None
                else None
            )
            return

        self._last_tm_cache_key = cache_key
        self._tm_value, tm_val, tm_tooltip = self._compute_tm_display(
            seq_val, name_val, show_temp, settings
        )
        self.tm_text.value = tm_val
        self.tm_container.tooltip = tm_tooltip

        self.tm_text.color = (
            tm_colour(self._tm_value, scheme)
            if self._tm_value is not None
            else None
        )

        if show_temp != getattr(self, "_last_show_temp", None):
            self._last_show_temp = show_temp
            self.tm_container.visible = show_temp
            self.tm_divider.visible = show_temp
            row_body: list[ft.Control] = [
                self.checkbox_container,
                self.active_divider,
                self.name_scroll,
                self.divider,
                self.seq_scroll,
            ]
            if show_temp:
                row_body.extend([self.tm_divider, self.tm_container])
            self._row_content = ft.Row(
                row_body,
                spacing=0,
                vertical_alignment=ft.CrossAxisAlignment.START,
            )
            self._row_gesture_detector.content = self._row_content
            self.content = ft.Row(
                [self.drag_handle, self._row_gesture_detector],
                spacing=0,
                vertical_alignment=ft.CrossAxisAlignment.START,
            )
