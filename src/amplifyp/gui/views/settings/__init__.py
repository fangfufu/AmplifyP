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

"""Settings sub-controls package."""

from amplifyp.gui.views.settings.appearance_tile import AppearanceTile
from amplifyp.gui.views.settings.base_score_tile import (
    BaseScoreTile,
    ScoreTable,
)
from amplifyp.gui.views.settings.dimer_tile import DimerTile
from amplifyp.gui.views.settings.replication_tile import ReplicationTile
from amplifyp.gui.views.settings.settings_view import SettingsView
from amplifyp.gui.views.settings.tm_tile import TmTile

__all__ = [
    "AppearanceTile",
    "BaseScoreTile",
    "DimerTile",
    "ReplicationTile",
    "ScoreTable",
    "SettingsView",
    "TmTile",
]
