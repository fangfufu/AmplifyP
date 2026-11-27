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
import contextvars
import importlib
import logging
import sys
import threading
from collections.abc import Awaitable, Callable
from pathlib import Path
from typing import Any, cast

import flet as ft

from amplifyp.version import (
    LATEST_RELEASE_API_URL,
    _get_sha,
    _parse_to_tuple,
    fetch_latest_release_version,
    get_full_sha,
    get_git_sha,
    is_newer_version,
    should_check_for_updates,
)
from amplifyp.version import (
    fetch_latest_release_version as _core_fetch_latest_release_version,
)
from amplifyp.version import (
    get_version as _core_get_version,
)

logger = logging.getLogger(__name__)


def is_pyodide_runtime() -> bool:
    """Return True when running under Pyodide / WebAssembly (Emscripten)."""
    return sys.platform == "emscripten" or "pyodide" in sys.modules


def can_use_background_threads() -> bool:
    """Return True when the runtime supports spawning OS threads.

    Pyodide / WebAssembly environments run single-threaded and raise
    ``RuntimeError("can't start new thread")`` if a thread is started.
    """
    return not is_pyodide_runtime()


async def run_in_background_thread[**P, T](
    func: Callable[P, T],
    *args: P.args,
    **kwargs: P.kwargs,
) -> T:
    """Run ``func`` in the default executor when threads are available.

    Falls back to calling ``func(*args, **kwargs)`` directly on
    single-threaded runtimes (such as Pyodide / WebAssembly) or if starting
    an executor thread fails with ``RuntimeError("can't start new thread")``.
    """
    if can_use_background_threads():
        loop = asyncio.get_running_loop()
        ctx = contextvars.copy_context()
        claim_lock = threading.Lock()
        claimed = False
        done_event = asyncio.Event()
        result_box: list[T] = []
        exc_box: list[BaseException] = []

        def _signal_done() -> None:
            try:
                loop.call_soon_threadsafe(done_event.set)
            except RuntimeError:
                pass

        def _run() -> T:
            nonlocal claimed
            with claim_lock:
                if claimed:
                    return cast("T", None)
                claimed = True
            try:
                res = ctx.run(func, *args, **kwargs)
            except BaseException as err:
                exc_box.append(err)
                _signal_done()
                raise
            result_box.append(res)
            _signal_done()
            return res

        try:
            future = loop.run_in_executor(None, _run)
        except RuntimeError as ex:
            if "can't start new thread" not in str(ex).lower():
                raise
            # ThreadPoolExecutor.submit() enqueues the work item before
            # spawning a new thread; coordinate with any existing worker that
            # may already have dequeued _run so func never executes twice.
            with claim_lock:
                already_claimed = claimed
                claimed = True
            if already_claimed:
                await done_event.wait()
                if exc_box:
                    raise exc_box[0] from None
                return result_box[0]
            logger.info(
                "Background thread unavailable (%s); running callable "
                "directly on event loop",
                ex,
            )
        else:
            return await future
    return func(*args, **kwargs)


async def fetch_latest_release_version_async() -> str | None:
    """Fetch the latest release version tag without blocking the event loop.

    Uses ``run_in_background_thread`` on desktop/server runtimes. On Pyodide /
    WebAssembly runtimes, uses ``pyodide.http.pyfetch`` asynchronously so
    network I/O never blocks the Python worker event loop or spawns threads.
    """
    # Check identity against the core function so tests or callers patching
    # ``amplifyp.gui.utils.system.fetch_latest_release_version`` are honoured
    # rather than silently bypassed by ``pyodide.http.pyfetch``.
    if (
        not can_use_background_threads()
        and fetch_latest_release_version is _core_fetch_latest_release_version
    ):
        try:
            pyodide_http = importlib.import_module("pyodide.http")
            pyfetch = getattr(pyodide_http, "pyfetch", None)
        except ImportError:
            pyfetch = None

        if callable(pyfetch):
            pyfetch_fn = cast("Callable[[str], Awaitable[Any]]", pyfetch)
            try:
                async with asyncio.timeout(5.0):
                    response = await pyfetch_fn(LATEST_RELEASE_API_URL)
                    if getattr(response, "ok", True):
                        data = await response.json()
                        if isinstance(data, dict):
                            tag_name = data.get("tag_name")
                            if isinstance(tag_name, str):
                                stripped = tag_name.strip()
                                if stripped:
                                    return stripped
                    else:
                        logger.warning(
                            "Failed to check for updates: HTTP %s",
                            getattr(response, "status", "unknown"),
                        )
            except Exception as ex:
                logger.warning("Failed to check for updates: %s", ex)
            return None

    return await run_in_background_thread(fetch_latest_release_version)


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
    "can_use_background_threads",
    "capture_all_views_async",
    "capture_view_screenshot_async",
    "confirm_dismiss",
    "confirm_exit",
    "confirm_exit_async",
    "fetch_latest_release_version",
    "fetch_latest_release_version_async",
    "get_full_sha",
    "get_git_sha",
    "get_version",
    "is_newer_version",
    "is_pyodide_runtime",
    "on_window_event",
    "restore_state_and_auto_close_async",
    "run_in_background_thread",
    "should_check_for_updates",
]
