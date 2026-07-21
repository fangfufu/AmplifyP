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

"""Consolidated keyboard navigation and GUI UI helper utilities."""

from __future__ import annotations

import asyncio
import threading
from collections.abc import Callable
from typing import TYPE_CHECKING, Any

import flet as ft

if TYPE_CHECKING:
    pass

from amplifyp.gui.colours import GUIColours

# ==============================================================================
# UI Helper Functions & Classes (formerly ui.py)
# ==============================================================================


def show_error_dialog(page: ft.Page, title: str, message: str) -> None:
    """Show an error dialog popup."""

    def close_dlg(e: ft.Event[ft.Control]) -> None:
        """Close the error dialog and update the page."""
        dialog.open = False
        page.update()

    def on_dismiss(e: ft.Event[ft.Control]) -> None:
        """Remove the dialog from the page overlay when dismissed."""
        if dialog in page.overlay:
            page.overlay.remove(dialog)
            page.update()

    dialog = ft.AlertDialog(
        title=ft.Text(title, color=GUIColours.ERROR_RED),
        content=ft.Text(message),
        actions=[ft.TextButton("OK", on_click=close_dlg)],  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        actions_alignment=ft.MainAxisAlignment.END,
        on_dismiss=on_dismiss,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
    )
    page.overlay.append(dialog)
    dialog.open = True
    page.update()


class Debouncer:
    """A thread-based debounce helper for delaying UI actions."""

    def __init__(self, delay_seconds: float = 0.15) -> None:
        """Initialize the Debouncer."""
        self.delay_seconds = delay_seconds
        self._timer: threading.Timer | None = None

    def trigger(self, callback: Callable[[], None]) -> None:
        """Trigger the callback after the specified delay."""
        self.cancel()

        self._timer = threading.Timer(self.delay_seconds, callback)
        self._timer.daemon = True
        try:
            self._timer.start()
        except RuntimeError:
            self._timer = None
            callback()

    def cancel(self) -> None:
        """Cancel any pending callback execution."""
        if self._timer is not None:
            self._timer.cancel()
            self._timer = None


def initialise_score_fields(
    settings_map: dict[str, Any],
    prefix: str,
    row_headers: list[str],
    col_headers: list[str],
    on_change_handler: Callable[[ft.Event[ft.Control] | None], None],
    font_size: int,
) -> None:
    """Initialise a grid of text fields for a score table in settings_map."""
    for r_char in row_headers:
        for c_char in col_headers:
            key = f"{prefix}_{r_char}_{c_char}"
            settings_map[key] = ft.TextField(
                value="0",
                on_change=on_change_handler,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
                text_align=ft.TextAlign.CENTER,
                dense=True,
                width=38,
                height=36,
                content_padding=4,
                text_style=ft.TextStyle(
                    color=GUIColours.DIAGRAM_BLACK, size=font_size
                ),
            )


class NotificationHelper:
    """Helper class to manage user notifications and messages."""

    def __init__(self, page: ft.Page) -> None:
        """Initialize the NotificationHelper."""
        self.page = page
        self._snack_bar = ft.SnackBar(ft.Text(""), open=False)
        self.page.overlay.append(self._snack_bar)

    def show_message(self, message: str) -> None:
        """Show a message to the user via a SnackBar."""
        self._snack_bar.content = ft.Text(message)
        self._snack_bar.open = True
        self.page.update()


class BorderedCheckbox(ft.Container):  # type: ignore[misc]
    """A checkbox wrapped in a container with a border matching input fields."""

    def __init__(
        self,
        label: str,
        value: bool = False,
        on_change: Callable[[ft.Event[ft.Control] | None], None] | None = None,
    ) -> None:
        """Initialize the BorderedCheckbox."""
        self.checkbox = ft.Checkbox(
            label=label,
            value=value,
            on_change=on_change,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )
        super().__init__(
            content=self.checkbox,
            border=ft.Border.all(1, GUIColours.OUTLINE),
            border_radius=5,
            padding=ft.Padding(10, 0, 10, 0),
            height=48,
            alignment=ft.Alignment(-1, 0),
        )

    @property
    def value(self) -> bool:
        """Get the value of the inner checkbox."""
        return bool(self.checkbox.value)

    @value.setter
    def value(self, val: bool | str) -> None:
        """Set the value of the inner checkbox."""
        if isinstance(val, str):
            self.checkbox.value = val.lower() == "true"
        else:
            self.checkbox.value = bool(val)


async def focus_async(res: Any) -> None:
    """Await a coroutine returned from a focus call (e.g. control.focus())."""
    await res


# ==============================================================================
# Keyboard Navigation Helper (formerly keyboard.py)
# ==============================================================================


def handle_keyboard_event(controller: Any, e: ft.KeyboardEvent) -> None:
    """Handle global keyboard events for primer navigation."""
    if (
        not controller.input_view
        or controller.view_container.content != controller.input_view
    ):
        return

    focused = controller.input_view._currently_focused_control
    if not (
        focused
        and isinstance(focused.data, dict)
        and "idx" in focused.data
        and "field" in focused.data
    ):
        return

    if not isinstance(focused, ft.TextField):
        return

    idx = focused.data["idx"]
    field = focused.data["field"]

    from amplifyp.gui.views.input.primer.row import PrimerRow

    target_field: ft.TextField | None = None

    if e.key == "Tab":
        controls = controller.input_view.primer_input.primers_list.controls
        if field == "name":
            for row in controls:
                if isinstance(row, PrimerRow) and row.idx == idx:
                    target_field = row.seq_field
                    target_field.selection = ft.TextSelection(
                        base_offset=0, extent_offset=0
                    )
                    target_field.update()
                    break
            else:
                return
        elif field == "seq":
            next_row = None
            for row in controls:
                if isinstance(row, PrimerRow) and row.idx == idx + 1:
                    next_row = row
                    break
            if next_row:
                target_field = next_row.name_field
                target_field.selection = ft.TextSelection(
                    base_offset=0, extent_offset=0
                )
                target_field.update()
            else:
                return
        else:
            return

        controller.input_view._skip_seq_focus_reset = True
        res = target_field.focus()
        if asyncio.iscoroutine(res):
            controller.page.run_task(focus_async, res)
        return

    if e.key in ("Arrow Left", "Arrow Right"):
        if field == "name" and e.key == "Arrow Right":
            cursor_pos = focused.data.get("cursor_pos", 0)
            if cursor_pos != len(focused.value or ""):
                return
            controls = controller.input_view.primer_input.primers_list.controls
            for row in controls:
                if isinstance(row, PrimerRow) and row.idx == idx:
                    target_field = row.seq_field
                    break
            else:
                return
            target_field.selection = ft.TextSelection(
                base_offset=0, extent_offset=0
            )
            target_field.update()
        elif field == "seq" and e.key == "Arrow Left":
            cursor_pos = focused.data.get("cursor_pos", 0)
            if cursor_pos != 0:
                return
            target_field = focused
            for row in controller.input_view.primer_input.primers_list.controls:
                if isinstance(row, PrimerRow) and row.idx == idx:
                    target_field = row.name_field
                    break
            else:
                return
            name_len = len(target_field.value or "")
            target_field.selection = ft.TextSelection(
                base_offset=name_len, extent_offset=name_len
            )
            target_field.update()
        else:
            return

        controller.input_view._skip_seq_focus_reset = True
        res = target_field.focus()
        if asyncio.iscoroutine(res):
            controller.page.run_task(focus_async, res)
        return

    if e.key not in ("Arrow Up", "Arrow Down"):
        return

    if e.key == "Arrow Down":
        next_idx = idx + 1
    else:
        next_idx = idx - 1

    controls = controller.input_view.primer_input.primers_list.controls
    target_row = None

    for row in controls:
        if isinstance(row, PrimerRow) and row.idx == next_idx:
            target_row = row
            break

    if target_row:
        target_field = (
            target_row.name_field if field == "name" else target_row.seq_field
        )

    if target_field is not None:
        current_cursor = focused.data.get("cursor_pos", 0)
        target_text = target_field.value or ""
        new_cursor = min(current_cursor, len(target_text))

        target_field.selection = ft.TextSelection(
            base_offset=new_cursor, extent_offset=new_cursor
        )
        try:
            target_field.update()
        except RuntimeError:
            pass

        controller.input_view._skip_seq_focus_reset = True
        res = target_field.focus()
        if asyncio.iscoroutine(res):
            controller.page.run_task(focus_async, res)
