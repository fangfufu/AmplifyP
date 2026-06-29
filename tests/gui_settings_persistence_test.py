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

"""Tests for GUI settings persistence on Desktop and Web platforms."""

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, mock_open, patch

import flet as ft
import yaml

from amplifyp.gui.settings import GUISettings


def test_get_config_path() -> None:
    """Test that OS-specific config paths are resolved correctly."""
    settings = GUISettings()

    # Test Windows
    with (
        patch("sys.platform", "win32"),
        patch.dict(
            "os.environ",
            {"APPDATA": "C:\\Users\\test\\AppData\\Roaming"},
            clear=True,
        ),
    ):
        path = settings._get_config_path()
        assert "AppData" in str(path)
        assert "AmplifyP" in str(path)
        assert path.name == "settings.yaml"

    # Test macOS
    with (
        patch("sys.platform", "darwin"),
        patch.dict("os.environ", {"HOME": "/Users/test"}, clear=True),
    ):
        path = settings._get_config_path()
        assert (
            path.as_posix()
            == "/Users/test/Library/Application Support/AmplifyP/settings.yaml"
        )

    # Test Linux with XDG_CONFIG_HOME
    with (
        patch("sys.platform", "linux"),
        patch.dict(
            "os.environ", {"XDG_CONFIG_HOME": "/custom/config"}, clear=True
        ),
    ):
        path = settings._get_config_path()
        assert path.as_posix() == "/custom/config/amplifyp/settings.yaml"

    # Test Linux without XDG_CONFIG_HOME
    with (
        patch("sys.platform", "linux"),
        patch.dict("os.environ", {"HOME": "/home/test"}, clear=True),
    ):
        path = settings._get_config_path()
        assert path.as_posix() == "/home/test/.config/amplifyp/settings.yaml"


def test_desktop_save_to_local() -> None:
    """Test saving settings to local YAML file on desktop."""
    settings = GUISettings()
    settings["primability_cutoff"] = "0.75"

    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False

    mock_path = MagicMock(spec=Path)
    mock_path.parent = MagicMock()

    m_open = mock_open()
    with (
        patch.object(settings, "_get_config_path", return_value=mock_path),
        patch("builtins.open", m_open),
    ):
        settings.save_to_local(mock_page)

    mock_path.parent.mkdir.assert_called_once_with(parents=True, exist_ok=True)

    # Filter out external background file opens (e.g. xonsh history json)
    config_open_call = None
    for call in m_open.call_args_list:
        if call[0] and call[0][0] == mock_path:
            config_open_call = call
            break
    assert config_open_call is not None, (
        f"Config file path {mock_path} was not opened"
    )
    assert config_open_call[0][1] == "w"
    assert config_open_call[1].get("encoding") == "utf-8"

    # Check if YAML was dumped (contains primability_cutoff)
    written_data = "".join(
        call[0][0] for call in m_open.return_value.write.call_args_list
    )
    parsed = yaml.safe_load(written_data)
    assert parsed["primability_cutoff"] == "0.75"


def test_desktop_load_from_local() -> None:
    """Test loading settings from local YAML file on desktop."""
    settings = GUISettings()
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False

    mock_path = MagicMock(spec=Path)
    mock_path.exists.return_value = True

    yaml_content = "primability_cutoff: '0.65'\ndark_mode: true"
    m_open = mock_open(read_data=yaml_content)

    with (
        patch.object(settings, "_get_config_path", return_value=mock_path),
        patch("builtins.open", m_open),
    ):
        settings.load_from_local(mock_page)

    assert settings["primability_cutoff"] == "0.65"
    assert settings["dark_mode"] is True


def test_desktop_load_from_local_corrupted() -> None:
    """Test loading from a corrupted YAML file falls back to defaults."""
    settings = GUISettings()
    settings["primability_cutoff"] = "0.99"

    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False

    mock_path = MagicMock(spec=Path)
    mock_path.exists.return_value = True

    # Invalid YAML syntax
    m_open = mock_open(read_data="{invalid: yaml:")

    with (
        patch.object(settings, "_get_config_path", return_value=mock_path),
        patch("builtins.open", m_open),
    ):
        # Should not raise exception
        settings.load_from_local(mock_page)

    # Since it failed, the settings should still be whatever they were
    assert settings["primability_cutoff"] == "0.99"


def test_web_save_and_load_client_storage() -> None:
    """Test saving and loading settings on web using page.client_storage."""
    settings = GUISettings()
    settings["primability_cutoff"] = "0.85"
    settings["dark_mode"] = "system"

    mock_page = MagicMock()
    mock_page.web = True
    mock_page.client_storage = MagicMock()

    # Set up client storage mock
    mock_storage = {}

    def mock_set(key: str, value: Any) -> None:
        mock_storage[key] = value

    def mock_get(key: str) -> Any:
        return mock_storage.get(key)

    def mock_contains(key: str) -> bool:
        return key in mock_storage

    mock_page.client_storage.set = MagicMock(side_effect=mock_set)
    mock_page.client_storage.get = MagicMock(side_effect=mock_get)
    mock_page.client_storage.contains_key = MagicMock(side_effect=mock_contains)

    # Test Save
    settings.save_to_local(mock_page)

    # Assert key exists in client storage with prefix
    assert "amplifyp.settings.primability_cutoff" in mock_storage
    assert mock_storage["amplifyp.settings.primability_cutoff"] == "0.85"
    assert mock_storage["amplifyp.settings.dark_mode"] == "system"

    # Test Load
    new_settings = GUISettings()
    new_settings.load_from_local(mock_page)

    assert new_settings["primability_cutoff"] == "0.85"
    assert new_settings["dark_mode"] == "system"
