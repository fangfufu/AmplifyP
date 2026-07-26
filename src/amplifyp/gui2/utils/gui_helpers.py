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

from collections.abc import Callable
from typing import TYPE_CHECKING, Any

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import (
    QCheckBox,
    QHBoxLayout,
    QLineEdit,
    QMessageBox,
    QSizePolicy,
    QWidget,
)

from amplifyp.gui2.colours import GUIColours

if TYPE_CHECKING:
    pass


# ==============================================================================
# UI Helper Functions & Classes
# ==============================================================================


def show_error_dialog(parent: Any, title: str, message: str) -> None:
    """Show an error dialog popup."""
    QMessageBox.critical(parent, title, message)


class Debouncer:
    """A thread-based debounce helper for delaying UI actions."""

    def __init__(self, delay_ms: int = 150) -> None:
        """Initialize the Debouncer."""
        self.delay_ms = delay_ms
        self._timer: QTimer | None = None

    def trigger(self, callback: Callable[[], None]) -> None:
        """Trigger the callback after the specified delay."""
        self.cancel()
        self._timer = QTimer()
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(callback)
        self._timer.start(self.delay_ms)

    def cancel(self) -> None:
        """Cancel any pending callback execution."""
        if self._timer is not None:
            self._timer.stop()
            self._timer = None


def initialise_score_fields(
    settings_map: dict[str, Any],
    prefix: str,
    row_headers: list[str],
    col_headers: list[str],
    on_change_handler: Callable[[str], None],
    font_size: int = 12,
) -> dict[str, QLineEdit]:
    """Initialise a grid of text fields for a score table."""
    from PySide6.QtCore import Qt

    fields: dict[str, QLineEdit] = {}
    for r_char in row_headers:
        for c_char in col_headers:
            key = f"{prefix}_{r_char}_{c_char}"
            field = QLineEdit("0")
            field.setAlignment(Qt.AlignmentFlag.AlignCenter)  # type: ignore[arg-type]  # Qt.AlignCenter
            field.setFixedWidth(38)
            field.setFixedHeight(36)
            field.textChanged.connect(lambda t, k=key: on_change_handler(k))
            settings_map[key] = field
            fields[key] = field
    return fields


class NotificationHelper:
    """Helper class to manage user notifications and messages."""

    def __init__(self, parent: Any) -> None:
        """Initialize the NotificationHelper."""
        self.parent = parent
        self._toast_widget: QWidget | None = None
        self._toast_timer: QTimer | None = None

    def show_message(self, message: str) -> None:
        """Show a message to the user via a toast notification."""
        if self._toast_timer is not None:
            self._toast_timer.stop()

        if self._toast_widget is None:
            self._toast_widget = QWidget(self.parent)
            self._toast_widget.setObjectName("toast")
            self._toast_widget.setStyleSheet(
                """
                #toast {
                    background-color: #323232;
                    color: #ffffff;
                    border-radius: 4px;
                    padding: 8px 16px;
                }
                """
            )
            layout = QHBoxLayout(self._toast_widget)
            layout.setContentsMargins(8, 4, 8, 4)
            self._label = QLineEdit(message)
            self._label.setReadOnly(True)
            self._label.setStyleSheet(
                "background: transparent; color: #ffffff; border: none;"
            )
            layout.addWidget(self._label)

        self._label.setText(message)
        self._toast_widget.adjustSize()

        parent_rect = self.parent.geometry()
        x = parent_rect.center().x() - self._toast_widget.width() // 2
        y = parent_rect.bottom() - self._toast_widget.height() - 20
        self._toast_widget.move(x, y)
        self._toast_widget.show()

        self._toast_timer = QTimer()
        self._toast_timer.setSingleShot(True)
        self._toast_timer.timeout.connect(self._hide_toast)
        self._toast_timer.start(3000)

    def _hide_toast(self) -> None:
        """Hide the toast notification."""
        if self._toast_widget is not None:
            self._toast_widget.hide()


class BorderedCheckbox(QWidget):
    """A checkbox wrapped in a container with a border."""

    def __init__(
        self,
        label: str,
        value: bool = False,
        on_change: Callable[[bool], None] | None = None,
    ) -> None:
        """Initialize the BorderedCheckbox."""
        super().__init__()
        self.checkbox = QCheckBox(label)
        self.checkbox.setChecked(value)
        if on_change:
            self.checkbox.stateChanged.connect(lambda s: on_change(bool(s)))

        layout = QHBoxLayout(self)
        layout.setContentsMargins(10, 0, 10, 0)
        layout.addWidget(self.checkbox)
        self.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.setFixedHeight(48)

        self._update_style()

    def _update_style(self) -> None:
        """Apply border styling."""
        outline = GUIColours.OUTLINE
        self.setStyleSheet(
            f"QFrame {{ border: 1px solid {outline}; border-radius: 5px; }}"
        )

    @property
    def value(self) -> bool:
        """Get the value of the inner checkbox."""
        return self.checkbox.isChecked()

    @value.setter
    def value(self, val: bool) -> None:
        """Set the value of the inner checkbox."""
        self.checkbox.setChecked(bool(val))


# ==============================================================================
# Keyboard Navigation Helper
# ==============================================================================


def handle_keyboard_event(
    controller: Any,
    key: str,
    field_name: str | None = None,
    idx: int | None = None,
    cursor_pos: int = 0,
) -> None:
    """Handle global keyboard events for primer navigation."""
    if (
        not controller.input_view
        or controller.current_view != controller.input_view
    ):
        return

    if idx is None or field_name is None:
        return

    from amplifyp.gui2.views.input.primer.row import PrimerRow

    controls = controller.input_view.primer_input.primers_list.controls

    if key == "Tab":
        if field_name == "name":
            for row in controls:
                if isinstance(row, PrimerRow) and row.idx == idx:
                    row.seq_field.setFocus()
                    row.seq_field.selectAll()
                    break
        elif field_name == "seq":
            for _i, row in enumerate(controls):
                if isinstance(row, PrimerRow) and row.idx == idx + 1:
                    row.name_field.setFocus()
                    row.name_field.selectAll()
                    break
        return

    if key in ("Arrow Left", "Arrow Right"):
        if field_name == "name" and key == "Arrow Right":
            focused_val = getattr(
                controller.input_view.primer_input, "focused_value", ""
            )
            if cursor_pos != len(focused_val or ""):
                return
            for row in controls:
                if isinstance(row, PrimerRow) and row.idx == idx:
                    row.seq_field.setFocus()
                    row.seq_field.selectAll()
                    break
        elif field_name == "seq" and key == "Arrow Left":
            if cursor_pos != 0:
                return
            for row in controls:
                if isinstance(row, PrimerRow) and row.idx == idx:
                    row.name_field.setFocus()
                    break
        return

    if key not in ("Arrow Up", "Arrow Down"):
        return

    next_idx = idx + 1 if key == "Arrow Down" else idx - 1
    if next_idx < 0 or next_idx >= len(controls):
        return

    for row in controls:
        if isinstance(row, PrimerRow) and row.idx == next_idx:
            target_field = (
                row.name_field if field_name == "name" else row.seq_field
            )
            target_field.setFocus()
            break
