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

"""PCR simulation view with diagram and results."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame,
    QLabel,
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


class PCRView(QWidget):
    """PCR simulation view with diagram panel and result cards."""

    def __init__(
        self,
        input_data: GUIInput,
        settings: GUISettings,
    ) -> None:
        """Initialize the PCR view."""
        super().__init__()
        self.input_data = input_data
        self.settings = settings
        self._pcr_result: Any = None
        self._build_ui()

    def _build_ui(self) -> None:
        """Build the PCR view UI."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)

        # Diagram panel placeholder
        self.diagram_label = QLabel("PCR Diagram will appear here")
        self.diagram_label.setAlignment(
            Qt.AlignmentFlag.AlignCenter | Qt.AlignmentFlag.AlignVCenter
        )  # Qt.AlignCenter
        self.diagram_label.setStyleSheet(
            f"background-color: {GUIColours.DIAGRAM_BG}; "
            f"color: {GUIColours.DIAGRAM_BLACK}; "
            "font-size: 14px; padding: 40px;"
        )
        self.diagram_label.setMinimumHeight(200)

        # Results area
        self.results_scroll = QScrollArea()
        self.results_scroll.setWidgetResizable(True)
        self.results_widget = QWidget()
        self.results_layout = QVBoxLayout(self.results_widget)
        self.results_layout.setContentsMargins(0, 0, 0, 0)
        self.results_scroll.setWidget(self.results_widget)

        # Splitter
        splitter = QSplitter()
        splitter.setOrientation(Qt.Orientation.Horizontal)
        splitter.addWidget(self.diagram_label)
        splitter.addWidget(self.results_scroll)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)

        layout.addWidget(splitter)

    def run_pcr(self, keep_cards: bool = False) -> bool:
        """Run PCR simulation."""
        from amplifyp.pcr import run_pcr

        primers = self.input_data.get_active_primers()
        if len(primers) < 1:
            return False

        template = self.input_data.template
        if not template.strip():
            return False

        try:
            result = run_pcr(
                template=template,
                circular=self.input_data.template_circular,
                primers=primers,
            )
            self._pcr_result = result
            self._update_ui()
            return True
        except Exception as e:
            self.diagram_label.setText(f"PCR Error: {e}")
            return False

    def _update_ui(self) -> None:
        """Update UI with PCR results."""
        if self._pcr_result is None:
            return

        # Clear existing cards
        while self.results_layout.count():
            item = self.results_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        # Add result cards
        if hasattr(self._pcr_result, "amplicons"):
            for i, amplicon in enumerate(self._pcr_result.amplicons[:100]):
                card = QFrame()
                card.setFrameShape(QFrame.Shape.StyledPanel)
                card.setStyleSheet(
                    f"background-color: {GUIColours.SURFACE_VARIANT}; "
                    "border-radius: 4px; padding: 8px;"
                )
                card_layout = QVBoxLayout(card)

                title = QLabel(
                    f"Amplicon {i + 1}: {getattr(amplicon, 'length', '?')} bp"
                )
                title.setStyleSheet(
                    f"font-weight: bold; color: {GUIColours.TEXT_ON_SURFACE};"
                )
                card_layout.addWidget(title)

                if hasattr(amplicon, "forward_primer"):
                    fp = getattr(amplicon, "forward_primer", None)
                    if fp:
                        fp_label = QLabel(f"Fwd: {fp}")
                        fp_label.setStyleSheet(
                            f"color: {GUIColours.FWD_PRIMER};"
                        )
                        card_layout.addWidget(fp_label)

                if hasattr(amplicon, "reverse_primer"):
                    rp = getattr(amplicon, "reverse_primer", None)
                    if rp:
                        rp_label = QLabel(f"Rev: {rp}")
                        rp_label.setStyleSheet(
                            f"color: {GUIColours.REV_PRIMER};"
                        )
                        card_layout.addWidget(rp_label)

                self.results_layout.insertWidget(i, card)

    def open_all_cards(self) -> None:
        """Open all result cards (expand)."""
