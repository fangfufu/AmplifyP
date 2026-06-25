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

import flet as ft

from amplifyp.gui.controller import GUIController


def main(
    page: ft.Page, state_file: str | None = None, auto_close: bool = False
) -> None:
    """Main entry point for the Flet application.

    Args:
        page: The Flet page instance.
        state_file: Optional path to a YAML state file to load on startup.
        auto_close: If True, quit automatically after rendering is complete.
    """
    controller = GUIController(
        page, state_file=state_file, auto_close=auto_close
    )
    controller.initialise()
