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

"""Main Flet application entry point."""

import flet as ft

from amplifyp.gui import main as app_main


def main(page: ft.Page) -> None:
    """Flet entry point - delegates to amplifyp.gui."""
    app_main(page)


if __name__ == "__main__":  # pragma: no cover
    ft.run(main, upload_dir="uploads")
