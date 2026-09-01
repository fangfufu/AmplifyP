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

"""Tests for the logging configuration of the GUI application."""

import logging
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from amplifyp.gui.logger import (
    _get_log_dir,
    get_default_log_file_path,
    initialise_logging,
    reconfigure_logging,
    resolve_log_file_path,
)


@pytest.fixture(autouse=True)  # type: ignore[untyped-decorator]
def reset_logging_flag() -> None:
    """Reset the internal logging initialization flag before each test."""
    import amplifyp.gui.logger

    amplifyp.gui.logger._logging_initialised = False
    amplifyp.gui.logger._current_file_path = None
    amplifyp.gui.logger._current_is_web = False

    # Clean up our own handlers from root logger to prevent cross-test pollution
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        if (
            type(handler) is logging.StreamHandler
            or isinstance(handler, logging.handlers.RotatingFileHandler)
            or type(handler) is logging.FileHandler
        ):
            root_logger.removeHandler(handler)
            handler.close()


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


def test_get_default_log_file_path() -> None:
    """Test that the default log file path is resolved correctly."""
    default_path = get_default_log_file_path()
    assert "app.log" in default_path
    assert "amplifyp" in default_path.lower()


def test_resolve_log_file_path() -> None:
    """Test path resolution with '(Default)' and custom paths."""
    # Test with '(Default)'
    resolved = resolve_log_file_path("(Default)")
    assert resolved == get_default_log_file_path()

    # Test with custom path
    custom_path = "/custom/path/to/my.log"
    resolved = resolve_log_file_path(custom_path)
    assert resolved == custom_path


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
    """Test logging setup when running in desktop mode with custom path."""
    custom_log_file = str(tmp_path / "test_app.log")

    initialise_logging(is_web=False, log_file_path=custom_log_file)

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
    assert (tmp_path / "test_app.log").exists()


def test_initialise_logging_no_rotation(tmp_path: Path) -> None:
    """Test logging setup with rotation disabled uses FileHandler."""
    custom_log_file = str(tmp_path / "test_no_rotation.log")

    initialise_logging(
        is_web=False,
        log_file_path=custom_log_file,
        log_rotation_enabled=False,
    )

    root_logger = logging.getLogger()
    console_handlers = [
        h for h in root_logger.handlers if type(h) is logging.StreamHandler
    ]
    rotating_handlers = [
        h
        for h in root_logger.handlers
        if type(h) is logging.handlers.RotatingFileHandler
    ]
    file_handlers = [
        h for h in root_logger.handlers if type(h) is logging.FileHandler
    ]

    assert len(console_handlers) == 1
    assert len(rotating_handlers) == 0
    assert len(file_handlers) == 1
    assert (tmp_path / "test_no_rotation.log").exists()


def test_reconfigure_logging_rotation_toggle(tmp_path: Path) -> None:
    """Test toggling rotation at runtime switches handler type."""
    custom_log_file = str(tmp_path / "test_rotation_toggle.log")

    # Start with rotation enabled
    initialise_logging(
        is_web=False,
        log_file_enabled=True,
        log_file_path=custom_log_file,
        log_rotation_enabled=True,
    )
    root_logger = logging.getLogger()
    rotating = [
        h
        for h in root_logger.handlers
        if type(h) is logging.handlers.RotatingFileHandler
    ]
    simple = [h for h in root_logger.handlers if type(h) is logging.FileHandler]
    assert len(rotating) == 1
    assert len(simple) == 0

    # Disable rotation
    reconfigure_logging(
        log_file_enabled=True,
        log_file_path=custom_log_file,
        log_rotation_enabled=False,
        is_web=False,
    )
    rotating = [
        h
        for h in root_logger.handlers
        if type(h) is logging.handlers.RotatingFileHandler
    ]
    simple = [h for h in root_logger.handlers if type(h) is logging.FileHandler]
    assert len(rotating) == 0
    assert len(simple) == 1

    # Re-enable rotation
    reconfigure_logging(
        log_file_enabled=True,
        log_file_path=custom_log_file,
        log_rotation_enabled=True,
        is_web=False,
    )
    rotating = [
        h
        for h in root_logger.handlers
        if type(h) is logging.handlers.RotatingFileHandler
    ]
    simple = [h for h in root_logger.handlers if type(h) is logging.FileHandler]
    assert len(rotating) == 1
    assert len(simple) == 0


def test_reconfigure_logging_max_bytes(tmp_path: Path) -> None:
    """Test that custom max bytes is applied to RotatingFileHandler."""
    custom_log_file = str(tmp_path / "test_max_bytes.log")
    custom_max = 1024 * 1024  # 1 MB

    initialise_logging(
        is_web=False,
        log_file_path=custom_log_file,
        log_rotation_max_bytes=custom_max,
    )

    root_logger = logging.getLogger()
    rotating = [
        h
        for h in root_logger.handlers
        if isinstance(h, logging.handlers.RotatingFileHandler)
    ]
    assert len(rotating) == 1
    assert rotating[0].maxBytes == custom_max


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


def test_initialise_logging_console_disabled() -> None:
    """Test logging setup with console output disabled."""
    initialise_logging(
        is_web=True,
        log_console_enabled=False,
    )

    root_logger = logging.getLogger()
    console_handlers = [
        h for h in root_logger.handlers if type(h) is logging.StreamHandler
    ]

    assert len(console_handlers) == 0


def test_reconfigure_logging_console_toggle() -> None:
    """Test toggling console output at runtime."""
    # Start with console enabled
    initialise_logging(is_web=True, log_console_enabled=True)
    root_logger = logging.getLogger()
    console_handlers = [
        h for h in root_logger.handlers if type(h) is logging.StreamHandler
    ]
    assert len(console_handlers) == 1

    # Disable console
    reconfigure_logging(is_web=True, log_console_enabled=False)
    console_handlers = [
        h for h in root_logger.handlers if type(h) is logging.StreamHandler
    ]
    assert len(console_handlers) == 0

    # Re-enable console
    reconfigure_logging(is_web=True, log_console_enabled=True)
    console_handlers = [
        h for h in root_logger.handlers if type(h) is logging.StreamHandler
    ]
    assert len(console_handlers) == 1


def test_reconfigure_logging_level_change() -> None:
    """Test changing log levels at runtime."""
    initialise_logging(is_web=True)

    # Check default levels
    amplifyp_logger = logging.getLogger("amplifyp.gui")
    flet_logger = logging.getLogger("flet")
    assert amplifyp_logger.level == logging.INFO
    assert flet_logger.level == logging.INFO

    # Change levels to valid values
    reconfigure_logging(
        log_level_amplifyp="WARNING",
        log_level_flet="WARNING",
        is_web=True,
    )
    assert amplifyp_logger.level == logging.WARNING
    assert flet_logger.level == logging.WARNING

    # Test that setting flet log level to DEBUG clamps to INFO
    reconfigure_logging(
        log_level_amplifyp="WARNING",
        log_level_flet="DEBUG",
        is_web=True,
    )
    assert flet_logger.level == logging.INFO


def test_reconfigure_logging_file_toggle(tmp_path: Path) -> None:
    """Test toggling file logging at runtime."""
    custom_log_file = str(tmp_path / "test_reconfigure.log")

    # Start with file logging enabled
    initialise_logging(
        is_web=False,
        log_file_enabled=True,
        log_file_path=custom_log_file,
    )
    root_logger = logging.getLogger()
    file_handlers = [
        h
        for h in root_logger.handlers
        if isinstance(h, logging.handlers.RotatingFileHandler)
    ]
    assert len(file_handlers) == 1

    # Disable file logging
    reconfigure_logging(
        log_file_enabled=False,
        is_web=False,
    )
    file_handlers = [
        h
        for h in root_logger.handlers
        if isinstance(h, logging.handlers.RotatingFileHandler)
    ]
    assert len(file_handlers) == 0

    # Re-enable file logging with new path
    new_path = str(tmp_path / "test_reconfigure2.log")
    reconfigure_logging(
        log_file_enabled=True,
        log_file_path=new_path,
        is_web=False,
    )
    file_handlers = [
        h
        for h in root_logger.handlers
        if isinstance(h, logging.handlers.RotatingFileHandler)
    ]
    assert len(file_handlers) == 1
    assert (tmp_path / "test_reconfigure2.log").exists()


def test_reconfigure_logging_invalid_level() -> None:
    """Test that invalid log level strings fall back to INFO."""
    initialise_logging(is_web=True)

    # Set invalid level
    reconfigure_logging(
        log_level_amplifyp="INVALID_LEVEL",
        is_web=True,
    )
    amplifyp_logger = logging.getLogger("amplifyp.gui")
    assert amplifyp_logger.level == logging.INFO


def test_initialise_logging_graceful_fallback() -> None:
    """Test that folder creation errors fall back to console gracefully."""
    # Mock Path.mkdir to raise PermissionError, simulating a read-only directory
    with patch.object(
        Path, "mkdir", side_effect=PermissionError("Permission denied")
    ):
        initialise_logging(
            is_web=False,
            log_file_path="/some/readonly/path/test.log",
        )

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


def test_initialise_logging_reapplies_settings() -> None:
    """Test that calling initialise_logging when already initialised
    reapplies settings.
    """
    initialise_logging(
        is_web=True,
        log_level_amplifyp="WARNING",
        log_level_flet="CRITICAL",
    )
    amplifyp_logger = logging.getLogger("amplifyp")
    amplifyp_gui_logger = logging.getLogger("amplifyp.gui")
    flet_logger = logging.getLogger("flet")
    assert amplifyp_logger.level == logging.WARNING
    assert amplifyp_gui_logger.level == logging.WARNING
    assert flet_logger.level == logging.CRITICAL

    # Re-initialise with different levels (flet log level should clamp to INFO)
    initialise_logging(
        is_web=True,
        log_level_amplifyp="ERROR",
        log_level_flet="DEBUG",
    )
    assert amplifyp_logger.level == logging.ERROR
    assert amplifyp_gui_logger.level == logging.ERROR
    assert flet_logger.level == logging.INFO


def test_logger_settings_path_and_apply_stored_settings() -> None:
    """Test _get_settings_path and _apply_stored_settings edge cases."""
    from unittest.mock import mock_open

    from amplifyp.gui.logger import (
        _apply_stored_settings,
        _get_log_dir,
        _get_settings_path,
    )

    # 1. _get_log_dir Windows without APPDATA
    with (
        patch("sys.platform", "win32"),
        patch.dict("os.environ", {}, clear=True),
    ):
        path = _get_log_dir()
        assert "Roaming" in str(path)

    # 2. _get_settings_path Windows with APPDATA
    with (
        patch("sys.platform", "win32"),
        patch.dict(
            "os.environ", {"APPDATA": r"C:\Users\test\AppData"}, clear=True
        ),
    ):
        path = _get_settings_path()
        assert str(path).startswith(r"C:\Users\test\AppData")

    # 3. _get_settings_path Windows without APPDATA
    with (
        patch("sys.platform", "win32"),
        patch.dict("os.environ", {}, clear=True),
    ):
        path = _get_settings_path()
        assert "Roaming" in str(path)

    # 4. _get_settings_path Darwin with and without HOME
    with (
        patch("sys.platform", "darwin"),
        patch.dict("os.environ", {"HOME": "/Users/tester"}, clear=True),
    ):
        path = _get_settings_path()
        assert path.as_posix() == (
            "/Users/tester/Library/Application Support/AmplifyP/settings.yaml"
        )

    with (
        patch("sys.platform", "darwin"),
        patch.dict("os.environ", {}, clear=True),
    ):
        path = _get_settings_path()
        assert "Application Support" in str(path)

    # 5. _get_settings_path Linux with XDG_CONFIG_HOME
    with (
        patch("sys.platform", "linux"),
        patch.dict(
            "os.environ", {"XDG_CONFIG_HOME": "/custom/xdg"}, clear=True
        ),
    ):
        path = _get_settings_path()
        assert path.as_posix() == "/custom/xdg/amplifyp/settings.yaml"

    # 6. _apply_stored_settings error loading and custom log path
    with (
        patch.object(Path, "exists", return_value=True),
        patch("builtins.open", side_effect=OSError("read error")),
    ):
        _apply_stored_settings()

    with (
        patch.object(Path, "exists", return_value=True),
        patch(
            "builtins.open",
            mock_open(read_data="log_file_path: /custom/log.txt"),
        ),
    ):
        _apply_stored_settings()
