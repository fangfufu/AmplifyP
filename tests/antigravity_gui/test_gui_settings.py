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

"""Test GUI settings using Amplify 4's settings module."""

import pytest

import amplifyp.antigravity_gui.gui as gui
from amplifyp.settings import GLOBAL_REPLICATION_SETTINGS, ReplicationSettings


def test_gui_settings_exists() -> None:
    """Verify that gui.settings exists and is of the correct type."""
    assert hasattr(gui, "settings"), (
        "gui module should have a 'settings' attribute"
    )
    assert isinstance(gui.settings, ReplicationSettings), (
        "gui.settings should be an instance of Settings"
    )


def test_gui_settings_initialization() -> None:
    """Verify that gui.settings is initialized with correct values."""
    # Direct comparison fails because helper classes don't implement __eq__ and
    # deepcopy creates new instances, so we check individual fields or values.
    # In the original script, we checked primability_cutoff and
    # stability_cutoff.
    assert gui.settings.primability_cutoff == pytest.approx(
        GLOBAL_REPLICATION_SETTINGS.primability_cutoff
    )
    assert gui.settings.stability_cutoff == pytest.approx(
        GLOBAL_REPLICATION_SETTINGS.stability_cutoff
    )


def test_gui_settings_modification() -> None:
    """Verify that gui.settings can be modified and persists."""
    gui.settings.primability_cutoff = 0.5
    assert gui.settings.primability_cutoff == pytest.approx(0.5)
