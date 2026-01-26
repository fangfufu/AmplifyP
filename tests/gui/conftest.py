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

"""Pytest configuration and fixtures for GUI tests."""

import copy
from collections.abc import Generator

import pytest

import amplifyp.gui as gui
from amplifyp.settings import GLOBAL_REPLICATION_SETTINGS


@pytest.fixture(autouse=True)  # type: ignore[untyped-decorator]
def reset_gui_settings() -> Generator[None, None, None]:
    """Fixture to reset gui.settings to defaults after each test.

    This ensures that any modifications to the global gui.settings object
    during a test do not leak into other tests.
    """
    yield
    gui.settings = copy.deepcopy(GLOBAL_REPLICATION_SETTINGS)
