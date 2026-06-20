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

_logging_initialised = False


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


def initialise_logging(is_web: bool = False) -> None:
    """Initialise logging configuration for the GUI app.

    Sets up a console handler and, if not running in a web context,
    a rotating file handler.

    Args:
        is_web: Whether the application is running in a web context.
    """
    global _logging_initialised
    if _logging_initialised:
        return

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    )
    console_handler.setFormatter(console_formatter)
    root_logger.addHandler(console_handler)

    # File handler (desktop mode only)
    if not is_web:
        try:
            log_dir = _get_log_dir()
            log_dir.mkdir(parents=True, exist_ok=True)
            log_file = log_dir / "app.log"

            file_handler = logging.handlers.RotatingFileHandler(
                log_file,
                maxBytes=5 * 1024 * 1024,  # 5 MB
                backupCount=3,
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
            # Fallback gracefully if directory creation or file writing fails
            # (e.g. permission issues in sandbox or restricted systems)
            logging.warning(
                "Could not initialise file logging: %s. "
                "Console logging will be used instead.",
                e,
            )

    _logging_initialised = True
