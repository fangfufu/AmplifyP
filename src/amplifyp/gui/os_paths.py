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

from amplifyp.os_paths import (
    get_app_config_dir,
    get_default_log_file_path,
    get_last_state_path,
    get_settings_yaml_path,
)

__all__ = [
    "get_app_config_dir",
    "get_default_log_file_path",
    "get_last_state_path",
    "get_settings_yaml_path",
]
