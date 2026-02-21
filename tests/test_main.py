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

"""Tests for the GUI app."""

from unittest.mock import MagicMock

import flet as ft

from main import main


def test_gui_main() -> None:
    """Test the main function of the GUI app."""
    mock_page = MagicMock(spec=ft.Page)

    main(mock_page)

    assert mock_page.title == "AmplifyP"
    assert mock_page.vertical_alignment == ft.MainAxisAlignment.START
    mock_page.add.assert_called()
