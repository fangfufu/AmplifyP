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

"""UI helper classes and widgets."""

import threading
from collections.abc import Callable
from typing import Any

import flet as ft


def show_error_dialog(page: ft.Page, title: str, message: str) -> None:
    """Show an error dialog popup.

    Creates and displays a modal AlertDialog with the given title and
    message, styled with the error colour.

    Args:
        page: The Flet page instance.
        title: The dialog title.
        message: The error message to display.
    """
    from amplifyp.gui.colours import GUIColours

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
        """Initialize the Debouncer.

        Args:
            delay_seconds: The delay in seconds before triggering the callback.
        """
        self.delay_seconds = delay_seconds
        self._timer: threading.Timer | None = None

    def trigger(self, callback: Callable[[], None]) -> None:
        """Trigger the callback after the specified delay.

        Any pending callback is cancelled.
        """
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
    """Initialise a grid of text fields for a score table in settings_map.

    Creates ft.TextField controls for each combination of row and column
    headers, storing them in settings_map with keys formatted as
    '{prefix}_{row}_{col}'.

    Args:
        settings_map: Dictionary to store the created TextField controls.
        prefix: The prefix for the field keys.
        row_headers: List of row header characters.
        col_headers: List of column header characters.
        on_change_handler: Callback function for field change events.
        font_size: The font size for the field text.
    """
    from amplifyp.gui.colours import GUIColours

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
    """Helper class to manage user notifications and messages.

    Wraps flet SnackBar usage to allow easy swapping to dialogues or other
    components.
    """

    def __init__(self, page: ft.Page) -> None:
        """Initialize the NotificationHelper.

        Args:
            page: The Flet page instance for displaying notifications.
        """
        self.page = page
        self._snack_bar = ft.SnackBar(ft.Text(""), open=False)
        self.page.overlay.append(self._snack_bar)

    def show_message(self, message: str) -> None:
        """Show a message to the user via a SnackBar.

        Updates the SnackBar content and opens it on the page overlay.

        Args:
            message: The message to display.
        """
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
        from amplifyp.gui.colours import GUIColours

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
