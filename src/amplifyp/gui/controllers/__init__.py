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

"""GUI controllers for theme, navigation, and auto-update management."""

from amplifyp.gui.controllers.navigation import NavigationManager
from amplifyp.gui.controllers.theme import ThemeManager
from amplifyp.gui.controllers.updater import UpdateManager

__all__ = ["NavigationManager", "ThemeManager", "UpdateManager"]
