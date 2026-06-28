# Copyright (C) 2026 Fufu Fang
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

"""Main Flet application logic."""

import logging

import flet as ft

from amplifyp.gui.controller import GUIController
from amplifyp.gui.logger import initialise_logging
from amplifyp.gui.settings import GUISettings

logger = logging.getLogger(__name__)


def main(
    page: ft.Page, state_file: str | None = None, auto_close: bool = False
) -> None:
    """Main entry point for the Flet application.

    Args:
        page: The Flet page instance.
        state_file: Optional path to a YAML state file to load on startup.
        auto_close: If True, quit automatically after rendering is complete.
    """
    gui_settings = GUISettings()
    gui_settings.load_from_local(page)

    initialise_logging(
        is_web=page.web,
        log_level_amplifyp=gui_settings.get("log_level_amplifyp", "DEBUG"),
        log_level_flet=gui_settings.get("log_level_flet", "INFO"),
        log_console_enabled=gui_settings.get("log_console_enabled", True),
        log_file_enabled=gui_settings.get("log_file_enabled", not page.web),
        log_file_path=gui_settings.get("log_file_path", "(Default)"),
        log_rotation_enabled=gui_settings.get("log_rotation_enabled", True),
        log_rotation_max_bytes=gui_settings.get(
            "log_rotation_max_bytes", 5242880
        ),
    )
    logger.info("Starting AmplifyP GUI application")
    controller = GUIController(
        page, state_file=state_file, auto_close=auto_close
    )
    controller.initialise()
