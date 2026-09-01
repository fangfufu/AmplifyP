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

"""Main Flet application logic."""

import logging
from typing import Any

import flet as ft

from amplifyp.gui.controller import GUIController
from amplifyp.gui.logger import initialise_logging
from amplifyp.gui.settings import GUISettings

logger = logging.getLogger(__name__)


def _patch_flet_session() -> None:
    """Patch Flet Session to safely handle invoke method results.

    Prevents crashing the message receive loop when invoke method results
    arrive for controls that were already unregistered or deleted.
    """
    try:
        from flet.messaging.session import (  # pyright: ignore[reportMissingTypeStubs]
            Session,
        )

        if getattr(Session, "_amplifyp_patched", False):
            return

        def safe_handle_invoke_method_results(
            self: Any,
            control_id: int,
            call_id: str,
            result: Any,
            error: str | None,
        ) -> None:
            """Safely handle invoke method results without crashing."""
            index = getattr(self, f"_{Session.__name__}__index", None)
            method_calls = getattr(
                self, f"_{Session.__name__}__method_calls", None
            )
            method_call_results = getattr(
                self, f"_{Session.__name__}__method_call_results", None
            )
            if index is not None and control_id in index:
                if method_calls is not None:
                    evt = method_calls.pop(call_id, None)
                    if evt is None:
                        return
                    if method_call_results is not None:
                        method_call_results[evt] = (result, error)
                    evt.set()
            else:
                if method_calls is not None:
                    evt = method_calls.pop(call_id, None)
                    if evt is not None:
                        if method_call_results is not None:
                            method_call_results[evt] = (result, error)
                        evt.set()
                logger.debug(
                    "Ignored invoke method result for unregistered control %s",
                    control_id,
                )

        Session.handle_invoke_method_results = safe_handle_invoke_method_results  # pyright: ignore[reportAttributeAccessIssue]
        Session._amplifyp_patched = True  # pyright: ignore[reportAttributeAccessIssue]
    except Exception as e:
        logger.debug("Failed to patch Flet session: %s", e)


_patch_flet_session()


def main(
    page: ft.Page,
    state_file: str | None = None,
    auto_close: bool = False,
    export_screenshots: bool = False,
    screenshots_dir: str | None = None,
    window_width: int | None = None,
    window_height: int | None = None,
) -> None:
    """Main entry point for the Flet application.

    Args:
        page: The Flet page instance.
        state_file: Optional path to a YAML state file to load on startup.
        auto_close: If True, quit automatically after rendering is complete.
        export_screenshots: Save PNG screenshots of views.
        screenshots_dir: Optional target directory for saved PNG screenshots.
        window_width: Optional application window width in pixels.
        window_height: Optional application window height in pixels.
    """
    gui_settings = GUISettings()
    gui_settings.load_from_local(page)

    initialise_logging(
        is_web=page.web,
        log_level_amplifyp=gui_settings.get("log_level_amplifyp", "DEBUG"),
        log_level_flet=gui_settings.get("log_level_flet", "INFO"),
        log_console_enabled=gui_settings.get("log_console_enabled", True),
        log_file_enabled=gui_settings.get("log_file_enabled", not page.web),
        log_file_path=gui_settings.get("log_file_path", "(Default)"),
        log_rotation_enabled=gui_settings.get("log_rotation_enabled", True),
        log_rotation_max_bytes=gui_settings.get(
            "log_rotation_max_bytes", 5242880
        ),
    )
    logger.info("Starting AmplifyP GUI application")
    controller = GUIController(
        page,
        state_file=state_file,
        auto_close=auto_close,
        export_screenshots=export_screenshots,
        screenshots_dir=screenshots_dir,
        window_width=window_width,
        window_height=window_height,
    )
    controller.initialise()
