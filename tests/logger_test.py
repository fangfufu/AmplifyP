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

"""Tests for the logging configuration of the GUI application."""

import logging
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from amplifyp.gui.logger import _get_log_dir, initialise_logging


@pytest.fixture(autouse=True)  # type: ignore[untyped-decorator]
def reset_logging_flag() -> None:
    """Reset the internal logging initialization flag before each test."""
    import amplifyp.gui.logger

    amplifyp.gui.logger._logging_initialised = False

    # Clean up our own handlers from root logger to prevent cross-test pollution
    root_logger = logging.getLogger()
    for handler in list(root_logger.handlers):
        if type(handler) is logging.StreamHandler or isinstance(
            handler, logging.handlers.RotatingFileHandler
        ):
            root_logger.removeHandler(handler)


def test_get_log_dir() -> None:
    """Test that OS-specific log directories are resolved correctly."""
    # Test Windows
    with (
        patch("sys.platform", "win32"),
        patch.dict(
            "os.environ",
            {"APPDATA": "C:\\Users\\test\\AppData\\Roaming"},
            clear=True,
        ),
    ):
        path = _get_log_dir()
        assert "AppData" in str(path)
        assert "AmplifyP" in str(path)
        assert path.name == "AmplifyP"

    # Test macOS
    with (
        patch("sys.platform", "darwin"),
        patch.dict("os.environ", {"HOME": "/Users/test"}, clear=True),
    ):
        path = _get_log_dir()
        assert (
            path.as_posix()
            == "/Users/test/Library/Application Support/AmplifyP"
        )

    # Test Linux with XDG_CONFIG_HOME
    with (
        patch("sys.platform", "linux"),
        patch.dict(
            "os.environ", {"XDG_CONFIG_HOME": "/custom/config"}, clear=True
        ),
    ):
        path = _get_log_dir()
        assert path.as_posix() == "/custom/config/amplifyp"

    # Test Linux without XDG_CONFIG_HOME
    with (
        patch("sys.platform", "linux"),
        patch.dict("os.environ", {"HOME": "/home/test"}, clear=True),
    ):
        path = _get_log_dir()
        assert path.as_posix() == "/home/test/.config/amplifyp"


def test_initialise_logging_web() -> None:
    """Test logging setup when running in web mode."""
    initialise_logging(is_web=True)

    root_logger = logging.getLogger()
    console_handlers = [
        h for h in root_logger.handlers if type(h) is logging.StreamHandler
    ]
    file_handlers = [
        h
        for h in root_logger.handlers
        if isinstance(h, logging.handlers.RotatingFileHandler)
    ]

    assert len(console_handlers) == 1
    assert len(file_handlers) == 0
    assert console_handlers[0].stream is sys.stdout


def test_initialise_logging_desktop(tmp_path: Path) -> None:
    """Test logging setup when running in desktop mode."""
    with (
        patch("amplifyp.gui.logger._get_log_dir", return_value=tmp_path),
    ):
        initialise_logging(is_web=False)

    root_logger = logging.getLogger()
    console_handlers = [
        h for h in root_logger.handlers if type(h) is logging.StreamHandler
    ]
    file_handlers = [
        h
        for h in root_logger.handlers
        if isinstance(h, logging.handlers.RotatingFileHandler)
    ]

    assert len(console_handlers) == 1
    assert len(file_handlers) == 1
    assert (tmp_path / "app.log").exists()


def test_initialise_logging_already_initialised() -> None:
    """Test that multiple calls to initialisation do not add handlers."""
    # First call
    initialise_logging(is_web=True)
    root_logger = logging.getLogger()
    console_handlers_1 = [
        h for h in root_logger.handlers if type(h) is logging.StreamHandler
    ]
    assert len(console_handlers_1) == 1

    # Second call
    initialise_logging(is_web=True)
    console_handlers_2 = [
        h for h in root_logger.handlers if type(h) is logging.StreamHandler
    ]
    assert len(console_handlers_2) == 1


def test_initialise_logging_graceful_fallback(tmp_path: Path) -> None:
    """Test that folder creation errors fall back to console gracefully."""
    # Mock mkdir to raise an exception (simulate permission error)
    mock_dir = MagicMock(spec=Path)
    mock_dir.mkdir.side_effect = PermissionError("Permission denied")

    with (
        patch("amplifyp.gui.logger._get_log_dir", return_value=mock_dir),
        patch("logging.warning") as mock_warning,
    ):
        initialise_logging(is_web=False)

    # Should fall back to console only
    root_logger = logging.getLogger()
    console_handlers = [
        h for h in root_logger.handlers if type(h) is logging.StreamHandler
    ]
    file_handlers = [
        h
        for h in root_logger.handlers
        if isinstance(h, logging.handlers.RotatingFileHandler)
    ]

    assert len(console_handlers) == 1
    assert len(file_handlers) == 0

    # Verify warning was logged
    mock_warning.assert_called_once()
