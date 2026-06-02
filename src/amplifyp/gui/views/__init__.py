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

"""GUI Views for the Flet application."""

from amplifyp.gui.state import GUIState

from .input_view import InputView
from .primer_dimer_view import PrimerDimerView
from .result_view import ResultView
from .settings_view import SettingsView

__all__ = [
    "GUIState",
    "InputView",
    "PrimerDimerView",
    "ResultView",
    "SettingsView",
]
