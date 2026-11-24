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

"""Tests for amplifyp.os_paths."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from amplifyp.os_paths import (
    get_app_config_dir,
    get_default_log_file_path,
    get_last_state_path,
    get_settings_yaml_path,
)


def test_windows_paths_with_appdata(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("APPDATA", r"C:\Users\Alice\AppData\Roaming")
    with patch("amplifyp.os_paths.sys.platform", "win32"):
        assert get_app_config_dir() == (
            Path(r"C:\Users\Alice\AppData\Roaming") / "AmplifyP"
        )
        assert get_settings_yaml_path() == (
            Path(r"C:\Users\Alice\AppData\Roaming")
            / "AmplifyP"
            / "settings.yaml"
        )
        assert get_default_log_file_path() == (
            Path(r"C:\Users\Alice\AppData\Roaming") / "AmplifyP" / "app.log"
        )
        assert get_last_state_path() == (
            Path(r"C:\Users\Alice\AppData\Roaming")
            / "AmplifyP"
            / "last_state.yaml"
        )


def test_windows_paths_without_appdata(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.delenv("APPDATA", raising=False)
    with (
        patch("amplifyp.os_paths.sys.platform", "win32"),
        patch(
            "amplifyp.os_paths.os.path.expanduser", return_value=str(tmp_path)
        ),
    ):
        assert get_app_config_dir() == (
            tmp_path / "AppData" / "Roaming" / "AmplifyP"
        )


def test_darwin_paths(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    with patch("amplifyp.os_paths.sys.platform", "darwin"):
        assert get_app_config_dir() == (
            tmp_path / "Library" / "Application Support" / "AmplifyP"
        )


def test_linux_paths_with_and_without_xdg(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    xdg = tmp_path / "xdg"
    monkeypatch.setenv("XDG_CONFIG_HOME", str(xdg))
    monkeypatch.setenv("HOME", str(tmp_path))
    with patch("amplifyp.os_paths.sys.platform", "linux"):
        assert get_app_config_dir() == xdg / "amplifyp"
        monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
        assert get_app_config_dir() == tmp_path / ".config" / "amplifyp"
