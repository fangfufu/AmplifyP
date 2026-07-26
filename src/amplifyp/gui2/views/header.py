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

"""Application header navigation control."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Any

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from amplifyp.gui2.colours import GUIColours
from amplifyp.gui2.settings import GUISettings
from amplifyp.gui2.utils.system import get_version

if TYPE_CHECKING:
    pass


class AppHeader(QWidget):
    """Application header component with navigation and state buttons."""

    def __init__(
        self,
        settings: GUISettings,
        on_switch_input: Callable[[], None],
        on_switch_settings: Callable[[], None],
        on_switch_about: Callable[[], None],
        on_pcr_click: Callable[[], None],
        on_dimers_click: Callable[[], None],
        on_save: Callable[[], bool],
        on_load: Callable[[], bool],
        pcr_button_ref: Any,
        dimers_button_ref: Any,
        on_clear_all: Callable[[], None] | None = None,
        on_switch_designer: Callable[[], None] | None = None,
        on_switch_designer_2d: Callable[[], None] | None = None,
    ) -> None:
        """Initialise the AppHeader navigation component."""
        super().__init__()
        self.settings = settings
        self.on_switch_input = on_switch_input
        self.on_switch_settings = on_switch_settings
        self.on_switch_about = on_switch_about
        self.on_pcr_click = on_pcr_click
        self.on_dimers_click = on_dimers_click
        self.on_save = on_save
        self.on_load = on_load
        self.on_clear_all = on_clear_all
        self.on_switch_designer = on_switch_designer
        self.on_switch_designer_2d = on_switch_designer_2d

        self.pcr_button_ref = pcr_button_ref
        self.dimers_button_ref = dimers_button_ref

        # Accept bool-returning callbacks by wrapping
        def _save_wrapper() -> None:
            on_save()

        def _load_wrapper() -> None:
            on_load()

        self._on_save = _save_wrapper
        self._on_load = _load_wrapper

        self._build_ui()

    def _build_ui(self) -> None:
        """Build the header UI."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 8, 16, 8)
        layout.setSpacing(8)

        # Title row
        title_layout = QHBoxLayout()
        title_layout.setSpacing(8)

        # Icon (if available)
        icon_path = self._get_icon_path()
        if icon_path:
            icon_label = QLabel()
            from PySide6.QtGui import QPixmap

            pixmap = QPixmap(icon_path)
            if not pixmap.isNull():
                icon_label.setPixmap(
                    pixmap.scaled(
                        32,
                        32,
                        aspectMode=getattr(Qt, "KeepAspectRatio", 0),  # type: ignore[arg-type]
                        transformationMode=getattr(
                            Qt, "SmoothTransformation", 1
                        ),  # type: ignore[arg-type]
                    )
                )
            title_layout.addWidget(icon_label)

        title_label = QLabel("AmplifyP")
        title_label.setStyleSheet(
            "font-size: 20px; font-weight: bold; color: "
            + GUIColours.TEXT_ON_SURFACE
        )
        title_layout.addWidget(title_label)

        title_layout.addSpacing(12)

        self.version_text = QLabel(get_version())
        self.version_text.setStyleSheet(
            f"font-size: 14px; color: {GUIColours.TEXT_ON_SURFACE}; "
            "opacity: 0.5;"
        )
        self.version_text.setObjectName("version_label")
        title_layout.addWidget(self.version_text)
        title_layout.addStretch()

        layout.addLayout(title_layout)

        # Button row
        btn_layout = QHBoxLayout()
        btn_layout.setSpacing(10)

        self._add_styled_button(btn_layout, "Input", self.on_switch_input)
        self.pcr_button_ref = self._add_styled_button(
            btn_layout, "PCR", self.on_pcr_click, enabled=False
        )
        self.dimers_button_ref = self._add_styled_button(
            btn_layout, "Primer Dimers", self.on_dimers_click, enabled=False
        )
        if self.on_switch_designer:
            self._add_styled_button(
                btn_layout, "Designer 1D", self.on_switch_designer
            )
        if self.on_switch_designer_2d:
            self._add_styled_button(
                btn_layout, "Designer 2D", self.on_switch_designer_2d
            )
        self._add_styled_button(btn_layout, "Settings", self.on_switch_settings)
        self._add_styled_button(btn_layout, "About", self.on_switch_about)

        # Divider
        divider = QFrame()
        divider.setFrameShape(QFrame.Shape.VLine)
        divider.setFrameShadow(QFrame.Shadow.Sunken)
        divider.setStyleSheet(f"color: {GUIColours.OUTLINE};")
        btn_layout.addWidget(divider)

        # Save/Clear/Load buttons
        self.save_btn_control = self._add_styled_button(
            btn_layout, "Save all", self._on_save
        )
        self.clear_btn_control = self._add_styled_button(
            btn_layout, "Clear all", self.on_clear_all or (lambda: None)
        )
        self.load_btn_control = self._add_styled_button(
            btn_layout, "Load all", self._on_load
        )

        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        # Initially hide save/clear/load
        self.save_btn_control.setVisible(False)
        self.clear_btn_control.setVisible(False)
        self.load_btn_control.setVisible(False)
        divider.setVisible(False)

        # Store divider reference for navigation manager
        self.header_divider = divider

    def _add_styled_button(
        self,
        layout: QHBoxLayout,
        text: str,
        callback: Callable[[], None],
        enabled: bool = True,
    ) -> QPushButton:
        """Add a styled button to the layout."""
        btn = QPushButton(text)
        btn.setEnabled(enabled)
        btn.setFixedHeight(32)
        btn.setStyleSheet(
            """
            QPushButton {
                background-color: #6750a4;
                color: #ffffff;
                border: none;
                border-radius: 4px;
                padding: 4px 12px;
                font-size: 13px;
            }
            QPushButton:hover {
                background-color: #7c67b5;
            }
            QPushButton:disabled {
                background-color: #9e9e9e;
                color: #bdbdbd;
            }
            """
        )
        btn.clicked.connect(callback)
        layout.addWidget(btn)
        return btn

    def _get_icon_path(self) -> str | None:
        """Get the icon path."""
        import os

        icon_path = os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "assets",
            "images",
            "favicon.png",
        )
        return icon_path if os.path.exists(icon_path) else None

    def set_update_available(self, new_version: str) -> None:
        """Update the version text to show that a new version is available."""
        current_version = get_version()
        self.version_text.setText(
            f"{current_version} (Update {new_version} available!)"
        )
        self.version_text.setStyleSheet(
            "font-size: 14px; color: "
            + GUIColours.UPDATE_AVAILABLE_COLOUR
            + "; opacity: 1.0;"
        )
        self.version_text.setToolTip("Click to open download page")
        self.version_text.mousePressEvent = lambda e: self._open_release_page()  # type: ignore[attr-defined]
        self.version_text.setCursor(  # type: ignore[union-attr]
            getattr(Qt, "PointingHandCursor", 12)  # type: ignore[arg-type]
        )

    def _open_release_page(self) -> None:
        """Open the release page in a browser."""
        import webbrowser

        webbrowser.open("https://github.com/fangfufu/AmplifyP/releases")
