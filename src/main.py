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
import os

import flet as ft

from amplifyp.gui import main as app_main

state_file: str | None = None
auto_close: bool = False


def main(page: ft.Page) -> None:
    """Flet entry point - delegates to amplifyp.gui."""
    app_main(page, state_file=state_file, auto_close=auto_close)


def cli(args_list: list[str] | None = None) -> None:
    """CLI entry point for argparse and running the Flet app."""
    global state_file, auto_close
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
        "--web",
        action="store_true",
        help="Launch in web browser mode",
    )
    parsed_args = parser.parse_args(args_list)
    if parsed_args.auto_close and not parsed_args.state:
        parser.error("--auto-close requires --state")
    state_file = parsed_args.state
    auto_close = parsed_args.auto_close

    assets_dir = os.path.join(os.path.dirname(__file__), "assets")
    view_mode = (
        ft.AppView.WEB_BROWSER if parsed_args.web else ft.AppView.FLET_APP
    )
    port_number = 34521 if parsed_args.web else 0
    ft.run(
        main,
        upload_dir="uploads",
        assets_dir=assets_dir,
        view=view_mode,
        port=port_number,
    )


if __name__ == "__main__":  # pragma: no cover
    cli()
