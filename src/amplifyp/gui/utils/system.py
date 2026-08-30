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
import json
import logging
import os
import re
import subprocess
import urllib.request
from importlib.metadata import PackageNotFoundError
from pathlib import Path
from typing import Any

import flet as ft

logger = logging.getLogger(__name__)


# ==============================================================================
# Git & Version Info (formerly git.py)
# ==============================================================================


def _get_sha(full: bool = False) -> str:
    """Retrieve the git commit SHA (either 40-char full or 7-char short)."""
    try:
        if full:
            from amplifyp.gui.git_sha import (  # pyright: ignore[reportMissingImports]
                GIT_FULL_SHA as imported_sha,  # pyright: ignore[reportMissingImports]
            )
        else:
            from amplifyp.gui.git_sha import (  # pyright: ignore[reportMissingImports]
                GIT_SHA as imported_sha,  # pyright: ignore[reportMissingImports]
            )

        if imported_sha and imported_sha != "unknown":
            return str(imported_sha)
    except ImportError:
        pass

    try:
        import js  # type: ignore[import-not-found, unused-ignore]

        if hasattr(js, "window") and hasattr(js.window, "__APP_SHA__"):
            sha = str(js.window.__APP_SHA__)
            if sha and sha != "unknown":
                return sha
    except (ImportError, AttributeError):
        pass

    try:
        cmd = (
            ["git", "rev-parse", "HEAD"]
            if full
            else ["git", "rev-parse", "--short", "HEAD"]
        )
        result = subprocess.run(  # noqa: S603
            cmd,
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass

    try:
        git_dir = os.path.join(
            os.path.dirname(
                os.path.dirname(
                    os.path.dirname(
                        os.path.dirname(
                            os.path.dirname(os.path.abspath(__file__))
                        )
                    )
                )
            ),
            ".git",
        )
        head_path = os.path.join(git_dir, "HEAD")
        if os.path.exists(head_path):
            with open(head_path, encoding="utf-8") as f:
                head_content = f.read().strip()
            if head_content.startswith("ref: refs/heads/"):
                ref_path = head_content.replace("ref: refs/heads/", "")
                ref_file = os.path.join(git_dir, ref_path)
                if os.path.exists(ref_file):
                    with open(ref_file, encoding="utf-8") as f:
                        full_sha = f.read().strip()
                    return full_sha if full else full_sha[:7]
            else:
                return head_content if full else head_content[:7]
    except OSError:
        pass

    try:
        dist_sha_path = os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "..", ".git-sha"
        )
        dist_sha_path = os.path.normpath(dist_sha_path)
        if os.path.exists(dist_sha_path):
            with open(dist_sha_path, encoding="utf-8") as f:
                return f.read().strip()
    except OSError:
        pass

    return "unknown"


def get_git_sha() -> str:
    """Get the short git commit SHA (7 chars), or 'unknown' if not available."""
    return _get_sha(full=False)


def get_full_sha() -> str:
    """Get the full git commit SHA (40 chars), or 'unknown' if not available."""
    return _get_sha(full=True)


def get_version() -> str:
    """Return version string like 'v0.0.1 (abc1234f)' or 'v0.0.1 (unknown)'."""
    local_logger = logging.getLogger(__name__)
    try:
        from amplifyp import __version__ as pkg_version
    except ImportError:
        try:
            from importlib.metadata import version

            pkg_version = version("amplifyp")
        except PackageNotFoundError:
            local_logger.debug("amplifyp package version not found")
            pkg_version = "unknown"

    git_sha = get_git_sha()
    return f"{pkg_version} ({git_sha})"


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


# ==============================================================================
# Version Check Helpers (formerly version_check.py)
# ==============================================================================


def fetch_latest_release_version() -> str | None:
    """Fetch the latest release version tag from GitHub API."""
    url = "https://api.github.com/repos/fangfufu/AmplifyP/releases/latest"
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "AmplifyP-Updater"},
    )
    try:
        with urllib.request.urlopen(req, timeout=5) as response:  # noqa: S310
            data = json.loads(response.read().decode("utf-8"))
            if isinstance(data, dict):
                tag_name = data.get("tag_name")
                if isinstance(tag_name, str):
                    return tag_name.strip()
    except Exception as e:
        logger.warning("Failed to check for updates: %s", e)
    return None


def _parse_to_tuple(v: str) -> tuple[int, ...]:
    """Parse a version string into a comparable tuple of integers."""
    if v.startswith("v"):
        v = v[1:]

    match = re.match(r"^([\d.]+)(.*)$", v)
    if not match:
        return (0, 0, 0, 0, 0, 0)

    release_str, pre_str = match.groups()
    release_list = [int(x) for x in release_str.split(".") if x.isdigit()]

    while release_list and release_list[-1] == 0:
        release_list.pop()

    while len(release_list) < 4:
        release_list.append(0)

    if not pre_str or pre_str.startswith("+"):
        pre_type = 4
        pre_num = 0
    else:
        pre_str = pre_str.lower()
        pre_word_match = re.search(r"([a-z]+)", pre_str)
        if pre_word_match:
            pre_word = pre_word_match.group(1)
            if pre_word in ("a", "alpha"):
                pre_type = 1
            elif pre_word in ("b", "beta"):
                pre_type = 2
            elif pre_word in ("rc", "preview"):
                pre_type = 3
            elif pre_word in ("post", "rev", "patch"):
                pre_type = 5
            else:
                pre_type = 0
        else:
            pre_type = 0

        num_match = re.search(r"(\d+)", pre_str)
        pre_num = int(num_match.group(1)) if num_match else 0

    return (*release_list, pre_type, pre_num)


def is_newer_version(latest_tag: str, current_ver: str) -> bool:
    """Compare latest tag with current version to see if update is available."""
    latest_tuple = _parse_to_tuple(latest_tag)
    current_tuple = _parse_to_tuple(current_ver)

    return latest_tuple > current_tuple


def should_check_for_updates(
    frequency: str, last_check: float, current_time: float
) -> bool:
    """Determine if update checking is due based on configured frequency."""
    if frequency == "At Startup":
        return True
    if frequency == "Disabled":
        return False

    elapsed = current_time - last_check

    if frequency == "Once per Day":
        return elapsed >= 86400.0
    if frequency == "Once per Week":
        return elapsed >= 604800.0
    if frequency == "Once per Month":
        return elapsed >= 2592000.0

    return False
