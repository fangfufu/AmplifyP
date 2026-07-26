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

"""Primer dimer analysis view."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from PySide6.QtWidgets import (
    QFrame,
    QLabel,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from amplifyp.gui2.colours import GUIColours
from amplifyp.gui2.settings import GUISettings
from amplifyp.gui2.user_data import GUIInput

if TYPE_CHECKING:
    pass


class DimerView(QWidget):
    """Primer dimer analysis view."""

    def __init__(
        self,
        input_data: GUIInput,
        settings: GUISettings,
    ) -> None:
        """Initialize the Dimer view."""
        super().__init__()
        self.input_data = input_data
        self.settings = settings
        self._dimer_result: Any = None
        self._build_ui()

    def _build_ui(self) -> None:
        """Build the dimer view UI."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)

        self.results_scroll = QScrollArea()
        self.results_scroll.setWidgetResizable(True)
        self.results_widget = QWidget()
        self.results_layout = QVBoxLayout(self.results_widget)
        self.results_layout.setContentsMargins(0, 0, 0, 0)
        self.results_scroll.setWidget(self.results_widget)

        layout.addWidget(self.results_scroll)

    def run_analysis(self) -> bool:
        """Run primer dimer analysis."""
        from amplifyp.dimer import PrimerDimerGenerator

        primers = self.input_data.get_active_primers()
        if len(primers) < 1:
            return False

        try:
            pd_settings = self.settings.get_primer_dimer_settings()
            generator = PrimerDimerGenerator(
                primers=primers,
                settings=pd_settings,
                max_results=self.settings.get("MAX_DIMERS_RENDER", 100),
            )
            self._dimer_result = generator.run()
            self._update_ui()
            return True
        except Exception:
            self._clear_results()
            return False

    def _update_ui(self) -> None:
        """Update UI with dimer results."""
        self._clear_results()

        if self._dimer_result is None:
            return

        if hasattr(self._dimer_result, "dimers"):
            for i, dimer in enumerate(self._dimer_result.dimers[:100]):
                card = QFrame()
                card.setFrameShape(QFrame.Shape.StyledPanel)
                card.setStyleSheet(
                    f"background-color: {GUIColours.SURFACE_VARIANT}; "
                    "border-radius: 4px; padding: 8px;"
                )
                card_layout = QVBoxLayout(card)

                title = QLabel(f"Dimer {i + 1}: {getattr(dimer, 'type', '?')}")
                title.setStyleSheet(
                    f"font-weight: bold; color: {GUIColours.TEXT_ON_SURFACE};"
                )
                card_layout.addWidget(title)

                if hasattr(dimer, "tm"):
                    tm_label = QLabel(f"Tm: {getattr(dimer, 'tm', '?'):.1f}°C")
                    tm_label.setStyleSheet(
                        f"color: {GUIColours.TEXT_ON_SURFACE};"
                    )
                    card_layout.addWidget(tm_label)

                self.results_layout.insertWidget(i, card)

    def _clear_results(self) -> None:
        """Clear result cards."""
        while self.results_layout.count():
            item = self.results_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
