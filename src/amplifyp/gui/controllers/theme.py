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

import flet as ft

from amplifyp.gui.colours import GUIColours


class ThemeManager:
    """Manages application appearance and theme transition handling."""

    def __init__(self, controller: Any) -> None:
        """Initialize ThemeManager with a reference to the main controller."""
        self.controller = controller

    def apply_theme(self) -> None:
        """Apply theme settings (light/dark/system mode) to the page."""
        dark_mode_setting = self.controller.settings.get("dark_mode", False)
        is_dark = False
        if str(dark_mode_setting).lower() == "system":
            self.controller.page.theme_mode = ft.ThemeMode.SYSTEM
            self.controller.page.bg_color = None  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            is_dark = (
                str(self.controller.page.platform_brightness).lower() == "dark"
            )
        elif bool(dark_mode_setting) and str(dark_mode_setting).lower() not in (
            "false",
            "0",
            "no",
        ):
            self.controller.page.theme_mode = ft.ThemeMode.DARK
            self.controller.page.bg_color = None  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            is_dark = True
        else:
            self.controller.page.theme_mode = ft.ThemeMode.LIGHT
            self.controller.page.bg_color = GUIColours.WHITE  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            is_dark = False
        GUIColours.dark_mode = is_dark
        if (
            hasattr(self.controller, "header_container")
            and self.controller.header_container
        ):
            self.controller.header_container.bgcolor = GUIColours.SURFACE

    def on_platform_brightness_change(
        self, _e: ft.ControlEvent | None = None
    ) -> None:
        """Handle system brightness shifts."""
        self.apply_theme()
        active_view = self.controller.view_container.content
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
        self.controller.page.update()
