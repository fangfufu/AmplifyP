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

"""OS-specific user data paths shared across the GUI modules.

This module depends only on the standard library so it can be imported
early (for example by :mod:`amplifyp.gui.logger`) without pulling in Flet
or the rest of the GUI package.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


def get_settings_yaml_path() -> Path:
    """Get the OS-specific path to the user settings.yaml file.

    Returns:
        Path object pointing to settings.yaml:
        - Windows: %APPDATA%/AmplifyP/settings.yaml
        - macOS: ~/Library/Application Support/AmplifyP/settings.yaml
        - Linux: $XDG_CONFIG_HOME/amplifyp/settings.yaml or
            ~/.config/amplifyp/settings.yaml
    """
    if sys.platform.startswith("win"):
        appdata = os.environ.get("APPDATA")
        if appdata:
            return Path(appdata) / "AmplifyP" / "settings.yaml"
        return (
            Path(os.path.expanduser("~"))
            / "AppData"
            / "Roaming"
            / "AmplifyP"
            / "settings.yaml"
        )
    elif sys.platform.startswith("darwin"):
        home = os.environ.get("HOME") or os.path.expanduser("~")
        return (
            Path(home)
            / "Library"
            / "Application Support"
            / "AmplifyP"
            / "settings.yaml"
        )
    else:
        xdg_config = os.environ.get("XDG_CONFIG_HOME")
        if xdg_config:
            return Path(xdg_config) / "amplifyp" / "settings.yaml"
        home = os.environ.get("HOME") or os.path.expanduser("~")
        return Path(home) / ".config" / "amplifyp" / "settings.yaml"
