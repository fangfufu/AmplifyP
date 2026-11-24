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

"""Consolidated system, version checking, and lifecycle utilities for GUI."""

from __future__ import annotations

import asyncio
import logging
from pathlib import Path
from typing import Any

import flet as ft

from amplifyp.version import (
    _get_sha,
    _parse_to_tuple,
    fetch_latest_release_version,
    get_full_sha,
    get_git_sha,
    is_newer_version,
    should_check_for_updates,
)
from amplifyp.version import (
    get_version as _core_get_version,
)

logger = logging.getLogger(__name__)


# ==============================================================================
# Git & Version Info (re-exported from amplifyp.version)
# ==============================================================================


def get_version() -> str:
    """Return version string like 'v0.0.1 (abc1234f)' or 'v0.0.1 (unknown)'."""
    return _core_get_version(git_sha_fn=get_git_sha)


# ==============================================================================
# Lifecycle & Exit Helpers (formerly lifecycle.py)
# ==============================================================================


def confirm_dismiss(controller: Any, _e: ft.ControlEvent) -> None:
    """Close close confirmation dialogue."""
    dialog = controller._confirm_dialog
    if dialog:
        dialog.open = False
        controller.page.update()


async def confirm_exit_async(controller: Any) -> None:
    """Asynchronously destroy the application window."""
    try:
        controller.save_last_state()
        await controller.page.window.destroy()
    except RuntimeError:
        logger.debug("Window already closed, skipping destroy")


def confirm_exit(controller: Any, _e: ft.ControlEvent) -> None:
    """Launch the async window destruction task."""
    controller.page.run_task(controller.confirm_exit_async)


async def capture_view_screenshot_async(
    page: ft.Page, view_name: str, screenshots_dir: Path | None = None
) -> Path:
    """Capture a PNG screenshot of the current page and save to directory.

    Args:
        page: The Flet page instance.
        view_name: The name identifier for the view (e.g. 'input_view').
        screenshots_dir: Target directory. Defaults to './screenshots'.

    Returns:
        Path object pointing to the written screenshot file.
    """
    if screenshots_dir is None:
        screenshots_dir = Path.cwd() / "screenshots"

    screenshots_dir.mkdir(parents=True, exist_ok=True)
    file_path = screenshots_dir / f"{view_name}.png"

    try:
        screenshot_bytes = await page.take_screenshot()
        if screenshot_bytes:
            file_path.write_bytes(screenshot_bytes)
            logger.info("Saved screenshot for %s to %s", view_name, file_path)
            return file_path
        logger.warning("Failed to capture screenshot bytes for %s", view_name)
        raise RuntimeError(
            f"Failed to capture screenshot bytes for {view_name}"
        )
    except Exception as ex:
        logger.exception("Error capturing screenshot for %s: %s", view_name, ex)
        raise


async def capture_all_views_async(controller: Any) -> None:
    """Switch across Input, PCR, and Dimer views and capture screenshots."""
    try:
        screenshots_dir = getattr(controller, "screenshots_dir", None)
        if not isinstance(screenshots_dir, (str, Path)):
            screenshots_dir = None
        await asyncio.sleep(1.0)
        controller.update_pcr_button_state(sync=True)

        # 1. Input View
        controller.switch_view(None, controller.input_view)
        controller.page.update()
        await asyncio.sleep(1.0)
        await capture_view_screenshot_async(
            controller.page, "input_view", screenshots_dir=screenshots_dir
        )

        # 2. PCR View
        controller.switch_view(None, controller.pcr_view)
        pcr_btn = controller.pcr_button_ref.current
        if pcr_btn and not pcr_btn.disabled:
            controller.pcr_view.run_pcr()
            controller.pcr_view.open_all_cards()
        controller.page.update()
        await asyncio.sleep(1.0)
        await capture_view_screenshot_async(
            controller.page, "pcr_view", screenshots_dir=screenshots_dir
        )

        # 3. Primer Dimer View
        dimers_btn = controller.dimers_button_ref.current
        if dimers_btn and not dimers_btn.disabled:
            controller.dimers_view.run_analysis()
        controller.switch_view(None, controller.dimers_view)
        controller.page.update()
        await asyncio.sleep(1.0)
        await capture_view_screenshot_async(
            controller.page,
            "primer_dimer_view",
            screenshots_dir=screenshots_dir,
        )

        if getattr(controller, "auto_close", False):
            await controller.confirm_exit_async()
    except Exception:
        logger.exception("Error during screenshot export sequence")
        if getattr(controller, "auto_close", False):
            try:
                await controller.confirm_exit_async()
            except RuntimeError:
                pass


async def restore_state_and_auto_close_async(controller: Any) -> None:
    """Restore state from file and run screenshot export or auto-close."""
    await asyncio.sleep(0)
    if controller.state_file:
        controller._restore_state_from_file(controller.state_file)
    if getattr(controller, "export_screenshots", False):
        await capture_all_views_async(controller)
    elif getattr(controller, "auto_close", False) and controller.state_file:
        await controller._auto_close_and_quit_delayed()


async def auto_close_and_quit_delayed(
    controller: Any, _event: ft.ControlEvent | None = None
) -> None:
    """Run PCR/dimer analysis then quit for performance regression testing."""
    try:
        await asyncio.sleep(0)

        controller.update_pcr_button_state(sync=False)

        pcr_btn = controller.pcr_button_ref.current
        if pcr_btn and not pcr_btn.disabled:
            controller.pcr_view.run_pcr()

        dimers_btn = controller.dimers_button_ref.current
        if dimers_btn and not dimers_btn.disabled:
            controller.dimers_view.run_analysis()

        controller.page.update()
        await asyncio.sleep(1)

        await controller.confirm_exit_async()
    except Exception:
        logger.exception("Error during auto-close sequence")
        try:
            await controller.confirm_exit_async()
        except RuntimeError:
            pass


def on_window_event(controller: Any, e: ft.WindowEvent) -> None:
    """Handle desktop window events, showing close confirmation dialog."""
    if (
        e.data == "close"
        or getattr(e, "type", None) == ft.WindowEventType.CLOSE
    ):
        dialog = controller._confirm_dialog
        msg = "Are you sure you want to close AmplifyP?"
        if not controller.settings.get("auto_reload_on_startup", True):
            msg += " Unsaved changes will be lost."

        if not dialog:
            dialog = ft.AlertDialog(
                modal=True,
                title=ft.Text("Confirm Exit"),
                content=ft.Text(msg),
                actions=[  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
                    ft.TextButton("Yes", on_click=controller.confirm_exit),  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
                    ft.TextButton("No", on_click=controller.confirm_dismiss),  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
                ],
                actions_alignment=ft.MainAxisAlignment.END,
            )
            controller._confirm_dialog = dialog
        else:
            dialog.content = ft.Text(msg)

        if dialog not in controller.page.overlay:
            controller.page.overlay.append(dialog)
        dialog.open = True
        controller.page.update()


__all__ = [
    "_get_sha",
    "_parse_to_tuple",
    "auto_close_and_quit_delayed",
    "capture_all_views_async",
    "capture_view_screenshot_async",
    "confirm_dismiss",
    "confirm_exit",
    "confirm_exit_async",
    "fetch_latest_release_version",
    "get_full_sha",
    "get_git_sha",
    "get_version",
    "is_newer_version",
    "on_window_event",
    "restore_state_and_auto_close_async",
    "should_check_for_updates",
]
