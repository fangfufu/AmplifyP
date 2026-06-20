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

"""Tests for the AboutView."""

from unittest.mock import MagicMock

import flet as ft

from amplifyp.gui.controller import GUIController
from amplifyp.gui.views.about.about_view import AboutView


def test_about_view_initialisation() -> None:
    """Test initialising the AboutView."""
    mock_page = MagicMock(spec=ft.Page)
    view = AboutView(mock_page)

    assert len(view.controls) > 0
    # Check that first control is the header Row
    assert isinstance(view.controls[0], ft.Row)


def test_controller_about_view_integration() -> None:
    """Test that GUIController initialises and routes to AboutView."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.window = MagicMock()
    mock_page.overlay = []

    controller = GUIController(mock_page)
    controller.initialise()

    assert controller.about_view is not None
    assert isinstance(controller.about_view, AboutView)

    # Verify transition works
    mock_event = MagicMock()
    controller.switch_view(mock_event, controller.about_view)
    assert controller.view_container.content == controller.about_view
