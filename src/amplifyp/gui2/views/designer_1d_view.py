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

"""1D primer designer view."""

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

from amplifyp.gui2.colours import GUIColours
from amplifyp.gui2.settings import GUISettings
from amplifyp.gui2.user_data import GUIInput

if TYPE_CHECKING:
    pass


class Designer1DView(QWidget):
    """1D primer designer view."""

    def __init__(
        self,
        input_data: GUIInput,
        settings: GUISettings,
    ) -> None:
        """Initialize the Designer 1D view."""
        super().__init__()
        self.input_data = input_data
        self.settings = settings
        self._designer_result: Any = None
        self._build_ui()

    def _build_ui(self) -> None:
        """Build the Designer 1D UI."""
        layout = QHBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)

        # Template input
        template_widget = QWidget()
        template_layout = QVBoxLayout(template_widget)
        template_layout.addWidget(QLabel("Template Sequence:"))
        self.template_edit = QPlainTextEdit()
        self.template_edit.setPlainText(self.input_data.template)
        template_layout.addWidget(self.template_edit)
        layout.addWidget(template_widget, 1)

        # Results
        results_widget = QWidget()
        results_layout = QVBoxLayout(results_widget)
        results_layout.addWidget(QLabel("Designed Primers:"))

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        self.results_area = QWidget()
        self.results_layout = QVBoxLayout(self.results_area)
        self.results_layout.setContentsMargins(0, 0, 0, 0)
        scroll.setWidget(self.results_area)
        results_layout.addWidget(scroll)

        run_btn = QPushButton("Run Designer 1D")
        run_btn.clicked.connect(self._run_designer)
        results_layout.addWidget(run_btn)
        layout.addWidget(results_widget, 2)

    def _run_designer(self) -> None:
        """Run the 1D primer designer."""
        template = self.template_edit.toPlainText()
        if not template.strip():
            return

        try:
            from amplifyp.primer_designer_1d import PrimerDesigner1D

            replication_settings = self.settings.get_replication_settings()
            designer = PrimerDesigner1D(
                template=template,
                replication_settings=replication_settings,
            )
            self._designer_result = designer.run()
            self._update_ui()
        except Exception as e:
            print(f"Designer 1D error: {e}")

    def _update_ui(self) -> None:
        """Update UI with designer results."""
        while self.results_layout.count():
            item = self.results_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        if self._designer_result is None:
            return

        if hasattr(self._designer_result, "primers"):
            for i, primer in enumerate(self._designer_result.primers[:100]):
                card = QFrame()
                card.setFrameShape(QFrame.Shape.StyledPanel)
                card.setStyleSheet(
                    f"background-color: {GUIColours.SURFACE_VARIANT}; "
                    "border-radius: 4px; padding: 8px;"
                )
                card_layout = QVBoxLayout(card)

                name = getattr(primer, "name", f"Primer {i + 1}")
                seq = getattr(primer, "sequence", "")
                quality = getattr(primer, "quality", 0)

                title = QLabel(f"{name}: {seq}")
                title.setStyleSheet(
                    f"font-weight: bold; color: {GUIColours.TEXT_ON_SURFACE};"
                )
                card_layout.addWidget(title)

                quality_label = QLabel(f"Quality: {round(quality)}")
                quality_label.setStyleSheet(
                    f"color: {GUIColours.TEXT_ON_SURFACE};"
                )
                card_layout.addWidget(quality_label)

                self.results_layout.insertWidget(i, card)
