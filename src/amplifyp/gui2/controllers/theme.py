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

"""Theme controller for managing dark/light modes and brightness shifts."""

from typing import Any

from amplifyp.gui2.colours import GUIColours


class ThemeManager:
    """Manages application appearance and theme transition handling."""

    def __init__(self, controller: Any) -> None:
        """Initialize ThemeManager with a reference to the main controller."""
        self.controller = controller

    def apply_theme(self) -> None:
        """Apply theme settings (light/dark/system mode)."""
        dark_mode_setting = self.controller.settings.get("dark_mode", False)
        is_dark = False
        if str(dark_mode_setting).lower() == "system":
            # System theme: detect from OS
            is_dark = self._detect_system_theme()
        elif str(dark_mode_setting).lower() not in (
            "false",
            "0",
            "no",
        ):
            is_dark = True
        else:
            is_dark = False

        GUIColours.dark_mode = is_dark

        if hasattr(self.controller, "header_container"):
            self.controller.header_container.setStyleSheet(
                f"background-color: {GUIColours.SURFACE};"
            )

    def _detect_system_theme(self) -> bool:
        """Detect if the system theme is dark."""
        try:
            from PySide6.QtGui import QPalette
            from PySide6.QtWidgets import QApplication

            app = QApplication.instance()
            if app is not None:  # type: ignore[possibly-undefined]
                palette = app.palette()  # type: ignore[union-attr]
                if palette is not None:
                    window_text = palette.color(QPalette.ColorRole.WindowText)
                    # Dark if window text is light (white-ish)
                    return window_text.lightness() > 128
        except Exception:  # noqa: S110
            pass
        return False

    def on_platform_brightness_change(self) -> None:
        """Handle system brightness shifts."""
        self.apply_theme()
        active_view = self.controller.current_view
        if active_view == self.controller.input_view:
            self.controller.input_view.update_ui()
        else:
            self.controller.input_view_dirty = True

        if active_view == self.controller.settings_view:
            self.controller.settings_view.update_ui()

        if active_view == self.controller.pcr_view:
            self.controller.pcr_view.run_pcr(keep_cards=True)
        elif active_view == self.controller.dimers_view:
            self.controller.dimers_view.run_analysis()
        if hasattr(self.controller, "main_window"):
            self.controller.main_window.update()
