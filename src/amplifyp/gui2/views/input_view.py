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

"""Input view for template and primer management."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Any

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from amplifyp.gui2.colours import GUIColours
from amplifyp.gui2.settings import GUISettings
from amplifyp.gui2.user_data import GUIInput

if TYPE_CHECKING:
    pass


class InputView(QWidget):
    """Main input view with template and primer panels."""

    def __init__(
        self,
        input_data: GUIInput,
        settings: GUISettings,
        on_change: Callable[[], None],
    ) -> None:
        """Initialize the Input view."""
        super().__init__()
        self.input_data = input_data
        self.settings = settings
        self.on_change = on_change
        self.validation_errors: list[dict[str, Any]] = []
        self.focused_primer_index: int | None = None
        self.selected_indices: set[int] = set()
        self._skip_seq_focus_reset = True

        # Primer input
        self.primer_input: Any = None

        self._build_ui()

    def _build_ui(self) -> None:
        """Build the input view UI."""
        layout = QHBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)

        # Template panel
        template_widget = QWidget()
        template_layout = QVBoxLayout(template_widget)

        color = GUIColours.TEXT_ON_SURFACE
        template_label = QLabel("DNA Template:")
        template_label.setStyleSheet(
            f"font-size: 14px; font-weight: bold; color: {color};"
        )
        template_layout.addWidget(template_label)

        self.template_edit = QPlainTextEdit()
        self.template_edit.setPlainText(self.input_data.template)
        self.template_edit.textChanged.connect(self._on_template_change)
        template_layout.addWidget(self.template_edit)

        # Circular checkbox
        self.circular_cb = QCheckBox("Circular template")
        self.circular_cb.setChecked(self.input_data.template_circular)
        self.circular_cb.stateChanged.connect(self._on_circular_change)
        template_layout.addWidget(self.circular_cb)

        template_layout.addStretch()

        # Primer panel
        primer_widget = QWidget()
        primer_layout = QVBoxLayout(primer_widget)

        primer_header = QHBoxLayout()
        primer_label = QLabel("Primers:")
        primer_label.setStyleSheet(
            f"font-size: 14px; font-weight: bold; color: {color};"
        )
        primer_header.addWidget(primer_label)
        primer_header.addStretch()
        primer_layout.addLayout(primer_header)

        # Primer list (scrollable)
        self.primer_list_widget = QWidget()
        self.primer_list_layout = QVBoxLayout(self.primer_list_widget)
        self.primer_list_layout.setContentsMargins(0, 0, 0, 0)
        self.primer_list_layout.setSpacing(4)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.primer_list_widget)
        primer_layout.addWidget(scroll)

        # Add primer button
        add_btn = QPushButton("+ Add Primer")
        add_btn.clicked.connect(self._add_primer)
        primer_layout.addWidget(add_btn)

        primer_layout.addStretch()

        # Splitter
        splitter = QSplitter()
        splitter.setOrientation(Qt.Orientation.Horizontal)  # Horizontal
        splitter.addWidget(template_widget)
        splitter.addWidget(primer_widget)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 2)

        layout.addWidget(splitter)

        # Add initial primer row
        self._add_primer()

    def _on_template_change(self) -> None:
        """Handle template text changes."""
        self.input_data.template = self.template_edit.toPlainText()
        self.on_change()

    def _on_circular_change(self, state: int) -> None:
        """Handle circular template checkbox change."""
        self.input_data.template_circular = state == 2
        self.on_change()

    def _add_primer(self) -> None:
        """Add a new primer row."""
        idx = len(self.input_data.primers)
        self.input_data.primers.append({"name": "", "seq": "", "active": False})

        row = self._create_primer_row(idx)
        self.primer_list_layout.addWidget(row)
        self.primer_list_layout.addStretch()

    def _create_primer_row(self, idx: int) -> QWidget:
        """Create a primer input row."""
        row_widget = QWidget()
        row_layout = QHBoxLayout(row_widget)
        row_layout.setContentsMargins(0, 0, 0, 0)
        row_layout.setSpacing(4)

        primer = self.input_data.primers[idx]

        # Active checkbox
        cb = QCheckBox()
        cb.setChecked(primer.get("active", False))
        cb.stateChanged.connect(lambda s, i=idx: self._on_primer_active(i, s))
        row_layout.addWidget(cb)

        # Name field
        name_edit = QLineEdit(primer.get("name", ""))
        name_edit.setMinimumWidth(100)
        name_edit.setPlaceholderText("Name")
        name_edit.textChanged.connect(
            lambda t, i=idx: self._on_primer_name(i, t)
        )
        row_layout.addWidget(name_edit)

        # Sequence field
        seq_edit = QLineEdit(primer.get("seq", ""))
        seq_edit.setMinimumWidth(200)
        seq_edit.setPlaceholderText("Sequence")
        seq_edit.setStyleSheet("font-family: monospace;")
        seq_edit.textChanged.connect(lambda t, i=idx: self._on_primer_seq(i, t))
        row_layout.addWidget(seq_edit)

        # Spacer
        row_layout.addStretch()

        # Delete button
        del_btn = QPushButton("X")
        del_btn.setFixedSize(24, 24)
        del_btn.clicked.connect(lambda _, i=idx: self._delete_primer(i))
        row_layout.addWidget(del_btn)

        return row_widget

    def _on_primer_name(self, idx: int, text: str) -> None:
        """Handle primer name change."""
        if idx < len(self.input_data.primers):
            self.input_data.primers[idx]["name"] = text
            self.on_change()

    def _on_primer_seq(self, idx: int, text: str) -> None:
        """Handle primer sequence change."""
        if idx < len(self.input_data.primers):
            self.input_data.primers[idx]["seq"] = text
            self.on_change()

    def _on_primer_active(self, idx: int, state: int) -> None:
        """Handle primer activation change."""
        if idx < len(self.input_data.primers):
            self.input_data.primers[idx]["active"] = state == 2
            self.on_change()

    def _delete_primer(self, idx: int) -> None:
        """Delete a primer row."""
        if len(self.input_data.primers) <= 1:
            return
        self.input_data.primers.pop(idx)
        self.update_ui()
        self.on_change()

    def sync_to_state(self) -> None:
        """Sync input data to UI controls before save."""
        self.template_edit.blockSignals(True)
        self.template_edit.setPlainText(self.input_data.template)
        self.template_edit.blockSignals(False)

        self.circular_cb.blockSignals(True)
        self.circular_cb.setChecked(self.input_data.template_circular)
        self.circular_cb.blockSignals(False)

        # Sync primer fields from data model to UI
        row_widgets = []
        for i in range(self.primer_list_layout.count()):
            item = self.primer_list_layout.itemAt(i)
            if item and item.widget() and i < len(self.input_data.primers):
                row_widgets.append((i, item.widget()))

        for idx, row_widget in row_widgets:
            primer = self.input_data.primers[idx]
            layout = row_widget.layout()
            if layout and layout.count() >= 3:
                # Col 1 = name field, Col 2 = seq field
                name_edit = layout.itemAt(1).widget()
                seq_edit = layout.itemAt(2).widget()
                if isinstance(name_edit, QLineEdit):
                    name_edit.blockSignals(True)
                    name_edit.setText(primer.get("name", ""))
                    name_edit.blockSignals(False)
                if isinstance(seq_edit, QLineEdit):
                    seq_edit.blockSignals(True)
                    seq_edit.setText(primer.get("seq", ""))
                    seq_edit.blockSignals(False)

    def update_ui(self) -> None:
        """Update UI from state."""
        self.template_edit.blockSignals(True)
        self.template_edit.setPlainText(self.input_data.template)
        self.template_edit.blockSignals(False)

        self.circular_cb.blockSignals(True)
        self.circular_cb.setChecked(self.input_data.template_circular)
        self.circular_cb.blockSignals(False)

    def get_primers(self) -> list[dict[str, Any]]:
        """Get current primer data."""
        return list(self.input_data.primers)
