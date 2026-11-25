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

"""OS-specific user configuration and data paths for AmplifyP."""

from __future__ import annotations

import os
import sys
from pathlib import Path

__all__ = [
    "get_app_config_dir",
    "get_default_log_file_path",
    "get_last_state_path",
    "get_settings_yaml_path",
]


def get_app_config_dir() -> Path:
    """Get the OS-specific directory for AmplifyP user configuration files.

    Returns:
        Path object pointing to the configuration directory:
        - Windows: ``%APPDATA%/AmplifyP``
        - macOS: ``~/Library/Application Support/AmplifyP``
        - Linux/POSIX: ``$XDG_CONFIG_HOME/amplifyp`` or ``~/.config/amplifyp``
    """
    if sys.platform.startswith("win"):
        appdata = os.environ.get("APPDATA")
        if appdata:
            return Path(appdata) / "AmplifyP"
        return (
            Path(os.path.expanduser("~")) / "AppData" / "Roaming" / "AmplifyP"
        )
    if sys.platform.startswith("darwin"):
        home = os.environ.get("HOME") or os.path.expanduser("~")
        return Path(home) / "Library" / "Application Support" / "AmplifyP"
    xdg_config = os.environ.get("XDG_CONFIG_HOME")
    if xdg_config:
        return Path(xdg_config) / "amplifyp"
    home = os.environ.get("HOME") or os.path.expanduser("~")
    return Path(home) / ".config" / "amplifyp"


def get_settings_yaml_path() -> Path:
    """Get the OS-specific path to the user ``settings.yaml`` file."""
    return get_app_config_dir() / "settings.yaml"


def get_default_log_file_path() -> Path:
    """Get the OS-specific default path to the application log file."""
    return get_app_config_dir() / "app.log"


def get_last_state_path() -> Path:
    """Get the OS-specific path to the saved ``last_state.yaml`` file."""
    return get_app_config_dir() / "last_state.yaml"
