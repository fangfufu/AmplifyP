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

"""2D primer designer view with heatmap grid."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from amplifyp.gui2.colours import GUIColours, designer_2d_colour
from amplifyp.gui2.settings import GUISettings
from amplifyp.gui2.user_data import GUIInput

if TYPE_CHECKING:
    pass


class Designer2DView(QWidget):
    """2D primer designer view with heatmap grid."""

    def __init__(
        self,
        input_data: GUIInput,
        settings: GUISettings,
    ) -> None:
        """Initialize the Designer 2D view."""
        super().__init__()
        self.input_data = input_data
        self.settings = settings
        self._designer_result: Any = None
        self._build_ui()

    def _build_ui(self) -> None:
        """Build the Designer 2D UI."""
        layout = QHBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)

        # Template input
        template_widget = QWidget()
        template_layout = QVBoxLayout(template_widget)
        template_layout.addWidget(QLabel("Forward Template:"))
        self.fwd_template = QPlainTextEdit()
        self.fwd_template.setPlainText(self.input_data.template)
        template_layout.addWidget(self.fwd_template)

        template_layout.addWidget(QLabel("Reverse Template:"))
        self.rev_template = QPlainTextEdit()
        template_layout.addWidget(self.rev_template)
        layout.addWidget(template_widget, 1)

        # Grid results
        grid_widget = QWidget()
        grid_layout = QVBoxLayout(grid_widget)
        grid_layout.addWidget(QLabel("Quality Grid:"))

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        self.grid_area = QWidget()
        self.grid_layout = QVBoxLayout(self.grid_area)
        self.grid_layout.setContentsMargins(0, 0, 0, 0)
        scroll.setWidget(self.grid_area)
        grid_layout.addWidget(scroll)

        run_btn = QPushButton("Run Designer 2D")
        run_btn.clicked.connect(self._run_designer)
        grid_layout.addWidget(run_btn)
        layout.addWidget(grid_widget, 2)

    def _run_designer(self) -> None:
        """Run the 2D primer designer."""
        fwd = self.fwd_template.toPlainText()
        if not fwd.strip():
            return

        try:
            from amplifyp.primer_designer_2d import PrimerDesigner2D

            replication_settings = self.settings.get_replication_settings()
            designer = PrimerDesigner2D(
                fwd_template=fwd,
                rev_template=self.rev_template.toPlainText(),
                replication_settings=replication_settings,
            )
            self._designer_result = designer.run()
            self._update_ui()
        except Exception as e:
            print(f"Designer 2D error: {e}")

    def _update_ui(self) -> None:
        """Update UI with designer 2D results."""
        while self.grid_layout.count():
            item = self.grid_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        if self._designer_result is None:
            return

        scheme = self.settings.get("designer_2d_colour_scheme", "Blue-Orange")

        if hasattr(self._designer_result, "grid"):
            grid = self._designer_result.grid
            if grid:
                min_val = (
                    min(min(row) for row in grid if row) if any(grid) else 0
                )
                max_val = (
                    max(max(row) for row in grid if row) if any(grid) else 1
                )

                for _row_idx, row in enumerate(grid):
                    row_widget = QHBoxLayout()
                    for _col_idx, value in enumerate(row):
                        cell = QFrame()
                        cell.setFixedSize(30, 30)
                        colour = designer_2d_colour(
                            value, min_val, max_val, scheme
                        )
                        if colour:
                            cell.setStyleSheet(
                                f"background-color: {colour}; "
                                "border: 1px solid " + GUIColours.OUTLINE + ";"
                            )
                        else:
                            bg = GUIColours.SURFACE_VARIANT
                            cell.setStyleSheet(
                                f"background-color: {bg}; "
                                "border: 1px solid " + GUIColours.OUTLINE + ";"
                            )
                        row_widget.addWidget(cell)
                    row_layout = QHBoxLayout()
                    row_layout.addLayout(row_widget)
                    self.grid_layout.addLayout(row_layout)
