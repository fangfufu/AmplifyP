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

"""GUI Views for the PySide6 application."""

from .about_view import AboutView
from .designer_1d_view import Designer1DView
from .designer_2d_view import Designer2DView
from .dimer_view import DimerView
from .input_view import InputView
from .pcr_view import PCRView
from .settings_view import SettingsView

__all__ = [
    "AboutView",
    "Designer1DView",
    "Designer2DView",
    "DimerView",
    "InputView",
    "PCRView",
    "SettingsView",
]
