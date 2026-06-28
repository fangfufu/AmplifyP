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

"""Logging configuration for the GUI application."""

import logging
import logging.handlers
import os
import sys
from pathlib import Path
from typing import Any

import yaml


def _get_log_dir() -> Path:
    """Get the OS-specific directory for log storage.

    Matches the behaviour in GUISettings._get_config_path.
    """
    if sys.platform.startswith("win"):
        appdata = os.environ.get("APPDATA")
        if appdata:
            return Path(appdata) / "AmplifyP"
        return (
            Path(os.path.expanduser("~")) / "AppData" / "Roaming" / "AmplifyP"
        )
    elif sys.platform.startswith("darwin"):
        home = os.environ.get("HOME") or os.path.expanduser("~")
        return Path(home) / "Library" / "Application Support" / "AmplifyP"
    else:
        xdg_config = os.environ.get("XDG_CONFIG_HOME")
        if xdg_config:
            return Path(xdg_config) / "amplifyp"
        home = os.environ.get("HOME") or os.path.expanduser("~")
        return Path(home) / ".config" / "amplifyp"


# Resolve default log file path at module import time (once, at startup)
_DEFAULT_LOG_FILE_PATH = str(_get_log_dir() / "app.log")

# Track current configuration state for reconfigure_logging()
_logging_initialised = False
_current_file_path: str | None = None
_current_is_web: bool = False


def _set_logger_levels(amplifyp_level: int, flet_level: int) -> None:
    """Set levels on all amplifyp and flet loggers."""
    if flet_level < logging.INFO:
        flet_level = logging.INFO
    for name in (
        "amplifyp",
        "amplifyp.gui",
        "src.amplifyp",
        "src.amplifyp.gui",
    ):
        logging.getLogger(name).setLevel(amplifyp_level)
    logging.getLogger("flet").setLevel(flet_level)
    for name in logging.root.manager.loggerDict:
        if name.startswith("flet"):
            logging.getLogger(name).setLevel(flet_level)


def _apply_stored_settings() -> None:
    """Apply stored logging settings from settings.yaml at import time.

    Reads the settings file directly to avoid circular imports,
    configuring logging so it works from the first module load.
    """
    global _logging_initialised, _current_file_path, _current_is_web

    # Read settings file directly to avoid circular imports
    settings_path = _get_settings_path()
    stored_settings: dict[str, Any] = {}
    if settings_path.exists():
        try:
            with open(settings_path, encoding="utf-8") as f:
                data = yaml.safe_load(f)
            if isinstance(data, dict):
                stored_settings = data
        except (OSError, yaml.YAMLError, ValueError):
            pass

    log_level_amplifyp = stored_settings.get("log_level_amplifyp", "INFO")
    log_level_flet = stored_settings.get("log_level_flet", "INFO")
    log_console_enabled = stored_settings.get("log_console_enabled", True)
    log_file_enabled = stored_settings.get("log_file_enabled", True)
    log_file_path = stored_settings.get("log_file_path", "(Default)")
    log_rotation_enabled = stored_settings.get("log_rotation_enabled", True)
    log_rotation_max_bytes = stored_settings.get(
        "log_rotation_max_bytes", 5242880
    )

    # Resolve path inline (avoid dependency on later-defined functions)
    if log_file_path == "(Default)":
        resolved_path = _DEFAULT_LOG_FILE_PATH
    else:
        resolved_path = log_file_path
    _current_file_path = resolved_path
    _current_is_web = False

    # Resolve levels inline
    amplifyp_level = getattr(logging, log_level_amplifyp.upper(), logging.INFO)
    flet_level = getattr(logging, log_level_flet.upper(), logging.INFO)

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    _set_logger_levels(amplifyp_level, flet_level)

    if log_console_enabled:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        )
        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)

    if log_file_enabled:
        try:
            log_dir = Path(resolved_path).parent
            log_dir.mkdir(parents=True, exist_ok=True)

            file_handler: logging.FileHandler
            if log_rotation_enabled:
                file_handler = logging.handlers.RotatingFileHandler(
                    resolved_path,
                    maxBytes=log_rotation_max_bytes,
                    backupCount=3,
                    encoding="utf-8",
                )
            else:
                file_handler = logging.FileHandler(
                    resolved_path,
                    encoding="utf-8",
                )
            file_handler.setLevel(logging.DEBUG)
            file_formatter = logging.Formatter(
                "%(asctime)s [%(levelname)s] %(name)s "
                "(%(filename)s:%(lineno)d): %(message)s"
            )
            file_handler.setFormatter(file_formatter)
            root_logger.addHandler(file_handler)
        except OSError:
            pass

    _logging_initialised = True


def _get_settings_path() -> Path:
    """Get the path to the settings.yaml file.

    Returns:
        Path object pointing to settings.yaml.
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


# Apply stored settings when the module is loaded
_apply_stored_settings()


def get_default_log_file_path() -> str:
    """Return the OS-specific default log file path.

    Returns:
        The default log file path string.
    """
    return _DEFAULT_LOG_FILE_PATH


def resolve_log_file_path(path: str) -> str:
    """Resolve a log file path, converting '(Default)' to the platform default.

    Args:
        path: The path string, possibly '(Default)'.

    Returns:
        The resolved absolute path string.
    """
    if path == "(Default)":
        return _DEFAULT_LOG_FILE_PATH
    return path


def _get_valid_level(level_str: str) -> int:
    """Get a valid logging level from a string, falling back to INFO.

    Args:
        level_str: The level string (e.g. 'DEBUG', 'INFO').

    Returns:
        The logging level constant.
    """
    level = getattr(logging, level_str.upper(), None)
    if level is not None and isinstance(level, int):
        return level
    return logging.INFO


def _find_console_handler(
    root_logger: logging.Logger,
) -> logging.StreamHandler[Any] | None:
    """Find the console (stdout) StreamHandler on the root logger.

    Args:
        root_logger: The root logger to search.

    Returns:
        The console StreamHandler, or None if not found.
    """
    for handler in root_logger.handlers:
        if isinstance(handler, logging.StreamHandler) and hasattr(
            handler, "stream"
        ):
            if handler.stream is sys.stdout:
                return handler
    return None


def _find_handler_by_type(
    root_logger: logging.Logger, handler_type: type
) -> logging.Handler | None:
    """Find a handler of a specific type on the root logger.

    Args:
        root_logger: The root logger to search.
        handler_type: The handler class to search for.

    Returns:
        The first matching handler, or None.
    """
    for handler in root_logger.handlers:
        if isinstance(handler, handler_type):
            return handler
    return None


def _remove_handlers_by_type(
    root_logger: logging.Logger, handler_type: type
) -> None:
    """Remove all handlers of a specific type from the root logger.

    Args:
        root_logger: The root logger to modify.
        handler_type: The handler class to remove.
    """
    root_logger.handlers = [
        h for h in root_logger.handlers if not isinstance(h, handler_type)
    ]


def initialise_logging(
    is_web: bool = False,
    log_level_amplifyp: str = "INFO",
    log_level_flet: str = "INFO",
    log_console_enabled: bool = True,
    log_file_enabled: bool = True,
    log_file_path: str = "(Default)",
    log_rotation_enabled: bool = True,
    log_rotation_max_bytes: int = 5242880,
) -> None:
    """Initialise logging configuration for the GUI app.

    Sets up a console handler and, if not running in a web context,
    a file handler (rotating or fixed depending on settings).

    Args:
        is_web: Whether the application is running in a web context.
        log_level_amplifyp: Log level for the amplifyp.gui logger.
        log_level_flet: Log level for the flet logger.
        log_console_enabled: Whether to enable console output.
        log_file_enabled: Whether to enable file logging.
        log_file_path: Path to the log file, or '(Default)' for the
            platform-specific default.
        log_rotation_enabled: Whether to use rotating file handler.
        log_rotation_max_bytes: Max bytes before rotation (only used
            if log_rotation_enabled is True).
    """
    global _logging_initialised, _current_file_path, _current_is_web

    if _logging_initialised:
        reconfigure_logging(
            log_level_amplifyp=log_level_amplifyp,
            log_level_flet=log_level_flet,
            log_console_enabled=log_console_enabled,
            log_file_enabled=log_file_enabled,
            log_file_path=log_file_path,
            is_web=is_web,
            log_rotation_enabled=log_rotation_enabled,
            log_rotation_max_bytes=log_rotation_max_bytes,
        )
        return

    _current_is_web = is_web
    _current_file_path = resolve_log_file_path(log_file_path)

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    _set_logger_levels(
        _get_valid_level(log_level_amplifyp),
        _get_valid_level(log_level_flet),
    )

    # Console handler
    if log_console_enabled:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        )
        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)

    # File handler (desktop mode only)
    if not is_web and log_file_enabled:
        try:
            log_dir = Path(_current_file_path).parent
            log_dir.mkdir(parents=True, exist_ok=True)

            file_handler: logging.FileHandler
            if log_rotation_enabled:
                file_handler = logging.handlers.RotatingFileHandler(
                    _current_file_path,
                    maxBytes=log_rotation_max_bytes,
                    backupCount=3,
                    encoding="utf-8",
                )
            else:
                file_handler = logging.FileHandler(
                    _current_file_path,
                    encoding="utf-8",
                )
            file_handler.setLevel(logging.DEBUG)
            file_formatter = logging.Formatter(
                "%(asctime)s [%(levelname)s] %(name)s "
                "(%(filename)s:%(lineno)d): %(message)s"
            )
            file_handler.setFormatter(file_formatter)
            root_logger.addHandler(file_handler)
        except OSError as e:
            logging.warning(
                "Could not initialise file logging: %s. "
                "Console logging will be used instead.",
                e,
            )

    _logging_initialised = True


def reconfigure_logging(
    log_level_amplifyp: str = "INFO",
    log_level_flet: str = "INFO",
    log_console_enabled: bool = True,
    log_file_enabled: bool = True,
    log_file_path: str = "(Default)",
    is_web: bool = False,
    log_rotation_enabled: bool = True,
    log_rotation_max_bytes: int = 5242880,
) -> None:
    """Reconfigure logging after initial setup.

    Called when user changes logging settings at runtime. Updates
    per-logger levels, adds/removes console and file handlers as needed.

    Args:
        log_level_amplifyp: Log level for the amplifyp.gui logger.
        log_level_flet: Log level for the flet logger.
        log_console_enabled: Whether to enable console output.
        log_file_enabled: Whether to enable file logging.
        log_file_path: Path to the log file, or '(Default)' for the
            platform-specific default.
        is_web: Whether the application is running in a web context.
        log_rotation_enabled: Whether to use rotating file handler.
        log_rotation_max_bytes: Max bytes before rotation (only used
            if log_rotation_enabled is True).
    """
    global _current_file_path, _current_is_web

    root_logger = logging.getLogger()

    # Update per-logger levels
    _set_logger_levels(
        _get_valid_level(log_level_amplifyp),
        _get_valid_level(log_level_flet),
    )

    # Handle console handler
    console_handler = _find_console_handler(root_logger)
    if log_console_enabled and console_handler is None:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        )
        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)
    elif not log_console_enabled and console_handler is not None:
        root_logger.removeHandler(console_handler)

    # Handle file handler
    resolved_path = resolve_log_file_path(log_file_path)
    old_path = _current_file_path
    _current_file_path = resolved_path
    _current_is_web = is_web

    # Remove existing file handlers if path changed, logging disabled,
    # web mode, or rotation type changed
    rotation_changed = False
    if old_path is not None and _current_file_path is not None:
        old_rotating = _find_handler_by_type(
            root_logger, logging.handlers.RotatingFileHandler
        )
        old_simple = _find_handler_by_type(root_logger, logging.FileHandler)
        if old_rotating is not None and not log_rotation_enabled:
            rotation_changed = True
        elif old_simple is not None and log_rotation_enabled:
            rotation_changed = True

    if old_path is not None and (
        resolved_path != old_path
        or not log_file_enabled
        or is_web
        or rotation_changed
    ):
        _remove_handlers_by_type(
            root_logger, logging.handlers.RotatingFileHandler
        )
        _remove_handlers_by_type(root_logger, logging.FileHandler)

    # Add file handler if enabled and not web
    if not is_web and log_file_enabled:
        existing_file = _find_handler_by_type(
            root_logger, logging.handlers.RotatingFileHandler
        )
        existing_simple = _find_handler_by_type(
            root_logger, logging.FileHandler
        )
        if existing_file is None and existing_simple is None:
            try:
                log_dir = Path(resolved_path).parent
                log_dir.mkdir(parents=True, exist_ok=True)

                file_handler: logging.FileHandler
                if log_rotation_enabled:
                    file_handler = logging.handlers.RotatingFileHandler(
                        resolved_path,
                        maxBytes=log_rotation_max_bytes,
                        backupCount=3,
                        encoding="utf-8",
                    )
                else:
                    file_handler = logging.FileHandler(
                        resolved_path,
                        encoding="utf-8",
                    )
                file_handler.setLevel(logging.DEBUG)
                file_formatter = logging.Formatter(
                    "%(asctime)s [%(levelname)s] %(name)s "
                    "(%(filename)s:%(lineno)d): %(message)s"
                )
                file_handler.setFormatter(file_formatter)
                root_logger.addHandler(file_handler)
            except OSError as e:
                logging.warning(
                    "Could not initialise file logging: %s. "
                    "Console logging will be used instead.",
                    e,
                )
