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
"""Main Flet application entry point."""

import sys
import traceback

import flet as ft

from amplifyp import main as _impl
from amplifyp.gui import main as app_main

state_file: str | None = None
auto_close: bool = False
export_screenshots: bool = False
screenshots_dir: str | None = None
window_width: int | None = None
window_height: int | None = None

__all__ = [
    "app_main",
    "auto_close",
    "cli",
    "export_screenshots",
    "main",
    "screenshots_dir",
    "state_file",
    "window_height",
    "window_width",
]


def main(page: ft.Page) -> None:
    """Flet entry point - delegates to amplifyp.gui."""
    _impl.state_file = state_file
    _impl.auto_close = auto_close
    _impl.export_screenshots = export_screenshots
    _impl.screenshots_dir = screenshots_dir
    _impl.window_width = window_width
    _impl.window_height = window_height
    _impl.main(page)


def cli(args_list: list[str] | None = None) -> None:
    """CLI entry point for argparse and running the Flet app."""
    global state_file, auto_close, export_screenshots, screenshots_dir
    global window_width, window_height
    _impl.cli(args_list)
    state_file = _impl.state_file
    auto_close = _impl.auto_close
    export_screenshots = _impl.export_screenshots
    screenshots_dir = _impl.screenshots_dir
    window_width = _impl.window_width
    window_height = _impl.window_height


if __name__ == "__main__":  # pragma: no cover
    try:
        cli()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
