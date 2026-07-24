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

import argparse
import logging
import os
import sys
import traceback

import flet as ft

from amplifyp.gui import main as app_main

state_file: str | None = None
auto_close: bool = False
export_screenshots: bool = False
screenshots_dir: str | None = None
window_width: int | None = None
window_height: int | None = None


def main(page: ft.Page) -> None:
    """Flet entry point - delegates to amplifyp.gui."""
    try:
        app_main(
            page,
            state_file=state_file,
            auto_close=auto_close,
            export_screenshots=export_screenshots,
            screenshots_dir=screenshots_dir,
            window_width=window_width,
            window_height=window_height,
        )
    except Exception:
        logging.getLogger(__name__).exception("Unhandled exception in main")
        raise


def cli(args_list: list[str] | None = None) -> None:
    """CLI entry point for argparse and running the Flet app."""
    global state_file, auto_close, export_screenshots, screenshots_dir
    global window_width, window_height
    parser = argparse.ArgumentParser(
        description="AmplifyP - Primer design and PCR simulation tool"
    )
    parser.add_argument(
        "-f",
        "--state",
        type=str,
        help="Path to a YAML state file to load on startup",
    )
    parser.add_argument(
        "--auto-close",
        action="store_true",
        help="Auto-quit after rendering completes (requires --state)",
    )
    parser.add_argument(
        "-s",
        "--screenshots",
        action="store_true",
        help="Save PNG screenshots of views (requires --state)",
    )
    parser.add_argument(
        "--screenshots-dir",
        type=str,
        help="Target directory for saved PNG screenshots",
    )
    parser.add_argument(
        "--window-width",
        type=int,
        help="Set application window width in pixels",
    )
    parser.add_argument(
        "--window-height",
        type=int,
        help="Set application window height in pixels",
    )
    parser.add_argument(
        "--web",
        action="store_true",
        help="Launch in web browser mode",
    )
    parsed_args = parser.parse_args(args_list)
    if parsed_args.auto_close and not parsed_args.state:
        parser.error("--auto-close requires --state")
    if parsed_args.screenshots and not parsed_args.state:
        parser.error("--screenshots requires --state")
    state_file = parsed_args.state
    auto_close = parsed_args.auto_close
    export_screenshots = parsed_args.screenshots
    screenshots_dir = parsed_args.screenshots_dir
    window_width = parsed_args.window_width
    window_height = parsed_args.window_height

    assets_dir = os.path.join(os.path.dirname(__file__), "assets")

    if sys.platform == "emscripten" or "pyodide" in sys.modules:
        view_mode = None
        port_number = 0
    else:
        view_mode = (
            ft.AppView.WEB_BROWSER if parsed_args.web else ft.AppView.FLET_APP
        )
        port_number = 34521 if parsed_args.web else 0

    ft.run(  # pyright: ignore[reportUnknownMemberType]
        main,
        upload_dir="uploads",
        assets_dir=assets_dir,
        view=view_mode,
        port=port_number,
    )


if __name__ == "__main__":  # pragma: no cover
    try:
        cli()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
