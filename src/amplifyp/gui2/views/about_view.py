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

"""About dialog for the PySide6 application."""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QTextBrowser,
    QVBoxLayout,
)

from amplifyp.gui2.colours import GUIColours
from amplifyp.gui2.settings import GUISettings
from amplifyp.gui2.utils.system import get_full_sha, get_git_sha, get_version

if TYPE_CHECKING:
    pass


class AboutView(QDialog):
    """About dialog showing application info."""

    def __init__(self, settings: GUISettings) -> None:
        """Initialize the About dialog."""
        super().__init__()
        self.settings = settings
        self.setWindowTitle("About AmplifyP")
        self.setMinimumSize(500, 400)
        self._build_ui()

    def _build_ui(self) -> None:
        """Build the about dialog UI."""
        layout = QVBoxLayout(self)
        layout.setSpacing(12)

        # Title
        color = GUIColours.TEXT_ON_SURFACE
        title = QLabel("AmplifyP")
        title.setStyleSheet(
            f"font-size: 24px; font-weight: bold; color: {color};"
        )
        layout.addWidget(title)

        # Version
        version_label = QLabel(get_version())
        version_label.setStyleSheet(
            f"font-size: 14px; color: {GUIColours.MUTED_GREY};"
        )
        layout.addWidget(version_label)

        # Description
        description = QTextBrowser()
        description.setPlainText(
            "AmplifyP is a modern primer design and PCR simulation tool.\n"
            "\n"
            "Features:\n"
            "- PCR simulation with amplicon prediction\n"
            "- Primer dimer analysis\n"
            "- 1D and 2D primer designer\n"
            "- Custom melting temperature calculation\n"
            "- Colour-blind friendly visualisation\n"
            "\n"
            "This is the PySide6 GUI (experimental).\n"
            "The original Flet-based GUI is available as the default."
        )
        description.setStyleSheet(
            f"font-size: 13px; color: {GUIColours.TEXT_ON_SURFACE};"
        )
        layout.addWidget(description)

        # Git info
        git_layout = QHBoxLayout()
        git_sha = QLabel(f"SHA: {get_git_sha()}")
        git_sha.setStyleSheet(
            f"font-size: 11px; color: {GUIColours.MUTED_GREY};"
        )
        full_sha = QLabel(f"Full SHA: {get_full_sha()}")
        full_sha.setStyleSheet(
            f"font-size: 11px; color: {GUIColours.MUTED_GREY};"
        )
        git_layout.addWidget(git_sha)
        git_layout.addWidget(full_sha)
        git_layout.addStretch()
        layout.addLayout(git_layout)

        # Close button
        close_btn = QPushButton("Close")
        close_btn.setFixedWidth(80)
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn, alignment=Qt.AlignmentFlag.AlignCenter)

    def _open_url(self, url: str) -> None:
        """Open a URL in the default browser."""
        import webbrowser

        webbrowser.open(url)
