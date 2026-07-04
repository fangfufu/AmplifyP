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

"""DiagnosticsTile component for Flet settings view."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.logger import get_default_log_file_path
from amplifyp.gui.settings import GUISettings

LOG_LEVELS = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]


class DiagnosticsTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Diagnostics settings."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.Event | None], None],
        header_size: int,
    ) -> None:
        """Initialise the DiagnosticsTile.

        Args:
            page: The Flet page instance.
            settings: The settings object.
            settings_map: A dictionary mapping setting keys to UI
                components for population and retrieval.
            on_change_handler: The handler to call when a setting changes.
            header_size: The size of the header text.
        """
        self._page = page
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler
        self.is_web = self._page.web

        from amplifyp.gui.utils.ui import BorderedCheckbox

        self.log_console_enabled = BorderedCheckbox(
            label="Console Output",
            value=settings.get("log_console_enabled", True),
            on_change=self.on_change_handler,
        )

        self.log_level_amplifyp = ft.Dropdown(
            label="AmplifyP Log Level",
            options=[ft.dropdown.Option(level) for level in LOG_LEVELS],
            value=settings.get("log_level_amplifyp", "INFO"),
            width=500,
            on_select=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        self.log_level_flet = ft.Dropdown(
            label="Flet Log Level",
            options=[
                ft.dropdown.Option(level)
                for level in LOG_LEVELS
                if level != "DEBUG"
            ],
            value=settings.get("log_level_flet", "INFO"),
            width=500,
            on_select=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        self.log_file_enabled = BorderedCheckbox(
            label="File Logging",
            value=settings.get("log_file_enabled", True),
            on_change=self._on_file_logging_change,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )

        self.log_rotation_enabled = BorderedCheckbox(
            label="Log Rotation",
            value=settings.get("log_rotation_enabled", True),
            on_change=self.on_change_handler,
        )

        current_max_mb = str(
            settings.get("log_rotation_max_bytes", 5242880) // (1024 * 1024)
        )
        mb_options = ["1", "5", "10", "20", "50"]
        if current_max_mb not in mb_options:
            mb_options.append(current_max_mb)
            mb_options.sort(key=int)

        self.log_rotation_max_bytes = ft.Dropdown(
            label="Max Log Size (MB)",
            options=[ft.dropdown.Option(mb) for mb in mb_options],
            value=current_max_mb,
            width=500,
            disabled=True,
            border_color=GUIColours.OUTLINE,
            on_select=self.on_change_handler,
        )

        current_path = settings.get("log_file_path", "(Default)")
        self.log_file_path = ft.Dropdown(
            label="Log File Path",
            width=500,
            on_select=self._on_log_file_path_change,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            border_color=GUIColours.OUTLINE,
        )
        self._update_log_file_path_options(current_path)

        self.web_info_text = ft.Text(
            "File logging is not available in web mode.",
            color=GUIColours.OUTLINE,
            italic=True,
            visible=self.is_web,
        )

        self.settings_map["log_console_enabled"] = self.log_console_enabled
        self.settings_map["log_level_amplifyp"] = self.log_level_amplifyp
        self.settings_map["log_level_flet"] = self.log_level_flet
        self.settings_map["log_file_enabled"] = self.log_file_enabled
        self.settings_map["log_file_path"] = self.log_file_path
        self.settings_map["log_rotation_enabled"] = self.log_rotation_enabled
        self.settings_map["log_rotation_max_bytes"] = (
            self.log_rotation_max_bytes
        )

        controls: list[ft.Control] = [
            self.log_console_enabled,
            self.log_level_amplifyp,
            self.log_level_flet,
            self.log_file_enabled,
            self.log_file_path,
            self.log_rotation_enabled,
            self.log_rotation_max_bytes,
        ]

        if self.is_web:
            controls.append(self.web_info_text)

        super().__init__(
            title=ft.Text(
                "Diagnostics Settings",
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Row(
                        [
                            ft.Container(
                                content=ft.Column(
                                    controls,
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=500,
                            ),
                        ],
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )

        self._update_file_path_field()
        self._update_rotation_controls()

    def _on_file_logging_change(self, e: ft.Event) -> None:
        """Handle file logging toggle - update rotation controls state."""
        self._update_rotation_controls()
        self.on_change_handler(e)

    def _update_rotation_controls(self) -> None:
        """Update rotation controls visibility/disabled state."""
        if self.is_web:
            self.log_rotation_enabled.visible = False
            self.log_rotation_max_bytes.visible = False
            return

        is_file_enabled = self.settings.get("log_file_enabled", True)
        self.log_rotation_enabled.visible = is_file_enabled
        self.log_rotation_max_bytes.visible = (
            is_file_enabled and self.settings.get("log_rotation_enabled", True)
        )
        self.log_rotation_max_bytes.disabled = (
            not is_file_enabled
            or not self.settings.get("log_rotation_enabled", True)
        )

    def _update_log_file_path_options(self, current_path: str) -> None:
        """Update dropdown options to include the current path and select it."""
        default_path = get_default_log_file_path()
        options = [
            ft.dropdown.Option(
                key="(Default)", text=f"(Default) {default_path}"
            ),
            ft.dropdown.Option(key="Select folder", text="Select folder"),
        ]
        if current_path != "(Default)" and current_path != "Select folder":
            options.insert(
                1, ft.dropdown.Option(key=current_path, text=current_path)
            )

        self.log_file_path.options = options
        self.log_file_path.value = current_path

    async def _on_log_file_path_change(self, e: ft.ControlEvent) -> None:
        if self.log_file_path.value == "Select folder":
            file_picker = ft.FilePicker()
            self._page.services.append(file_picker)
            self._page.update()
            try:
                selected_dir = await file_picker.get_directory_path(
                    dialog_title="Select Log Directory",
                )
                if selected_dir:
                    import os

                    new_path = os.path.join(selected_dir, "app.log")
                    self.settings["log_file_path"] = new_path
                    self._update_log_file_path_options(new_path)
                    self.on_change_handler(e)
                else:
                    prev_path = self.settings.get("log_file_path", "(Default)")
                    self._update_log_file_path_options(prev_path)
                    self._page.update()
            finally:
                if file_picker in self._page.services:
                    self._page.services.remove(file_picker)
                    self._page.update()
        else:
            self.on_change_handler(e)

    def _update_file_path_field(self) -> None:
        """Update the file path field based on current settings."""
        if self.is_web:
            self.log_file_path.visible = False
            self.log_file_enabled.visible = False
            return

        self.log_file_path.visible = True
        self.log_file_enabled.visible = True

        is_enabled = self.settings.get("log_file_enabled", True)
        self.log_file_path.disabled = not is_enabled

        current_path = self.settings.get("log_file_path", "(Default)")
        self._update_log_file_path_options(current_path)

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central settings."""
        self.log_console_enabled.value = bool(
            self.settings.get("log_console_enabled", True)
        )
        self.log_level_amplifyp.value = self.settings.get(
            "log_level_amplifyp", "DEBUG"
        )
        self.log_level_flet.value = self.settings.get("log_level_flet", "INFO")
        self.log_file_enabled.value = bool(
            self.settings.get("log_file_enabled", True)
        )
        self.log_rotation_enabled.value = bool(
            self.settings.get("log_rotation_enabled", True)
        )
        max_bytes = self.settings.get("log_rotation_max_bytes", 5242880)
        try:
            max_bytes = int(max_bytes)
        except (ValueError, TypeError):
            max_bytes = 5242880
        mb_value = str(max_bytes // (1024 * 1024))

        existing_keys = [opt.key for opt in self.log_rotation_max_bytes.options]
        if mb_value not in existing_keys:
            self.log_rotation_max_bytes.options.append(
                ft.dropdown.Option(mb_value)
            )
            self.log_rotation_max_bytes.options.sort(
                key=lambda opt: (
                    int(opt.key) if opt.key and opt.key.isdigit() else 0
                )
            )

        self.log_rotation_max_bytes.value = mb_value
        self._update_file_path_field()
        self._update_rotation_controls()
