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

"""GUI Views for the Flet application."""

from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput

from .about import AboutView
from .designer_1d import PrimerDesignerView
from .dimer import DimerView
from .input import InputView
from .pcr import PCRView
from .settings import SettingsView

__all__ = [
    "AboutView",
    "DimerView",
    "GUIInput",
    "GUISettings",
    "InputView",
    "PCRView",
    "PrimerDesignerView",
    "SettingsView",
]
