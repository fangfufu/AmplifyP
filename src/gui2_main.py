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
"""Main PySide6 application entry point."""

import argparse
import sys
import traceback

from amplifyp.gui2.app import main


def cli(args_list: list[str] | None = None) -> None:
    """CLI entry point for the PySide6 application."""
    parser = argparse.ArgumentParser(
        description=(
            "AmplifyP PySide6 GUI - Primer design and PCR simulation tool"
        )
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

    parsed_args = parser.parse_args(args_list)
    if parsed_args.auto_close and not parsed_args.state:
        parser.error("--auto-close requires --state")
    if parsed_args.screenshots and not parsed_args.state:
        parser.error("--screenshots requires --state")

    try:
        main(
            state_file=parsed_args.state,
            auto_close=parsed_args.auto_close,
            export_screenshots=parsed_args.screenshots,
            screenshots_dir=parsed_args.screenshots_dir,
            window_width=parsed_args.window_width,
            window_height=parsed_args.window_height,
        )
    except Exception:
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    cli()
