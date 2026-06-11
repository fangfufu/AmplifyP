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

"""Utility functions for sequence handling and state serialization."""

import os
import subprocess
import threading
from collections.abc import Callable
from typing import Any

import flet as ft
import yaml


def clean_sequence(seq: str) -> str:
    """Clean sequence of escaped and standard whitespaces."""
    if not seq:
        return ""
    clean = str(seq).replace("\\n", "").replace("\\t", "").replace("\\r", "")
    return "".join(clean.split())


def format_sequence(seq: str, wrap_length: int = 80) -> str:
    """Format sequence into lines of specified length."""
    clean = clean_sequence(seq)
    return "\n".join(
        [clean[i : i + wrap_length] for i in range(0, len(clean), wrap_length)]
    )


def serialize_state(state: dict[str, object]) -> str:
    """Serialize state dict to YAML string, handling multiline strings."""

    def multiline_presenter(dumper: yaml.Dumper, data: str) -> yaml.ScalarNode:
        if "\n" in data:
            return dumper.represent_scalar(
                "tag:yaml.org,2002:str", data, style="|"
            )
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    class _StateDumper(yaml.Dumper):
        pass

    _StateDumper.add_representer(str, multiline_presenter)
    return yaml.dump(state, Dumper=_StateDumper, sort_keys=False)


def create_overlapped_sequence_view(
    top_line: str,
    mid_line: str,
    bottom_line: str,
    font_family: str = "Roboto Mono",
    font_size: int = 14,
) -> ft.Text:
    """Create a Flet Text control showing visually aligned sequences.

    Uses TextSpans for the visual representation.
    """
    from amplifyp.gui.settings import GUIColors

    return ft.Text(
        spans=[
            ft.TextSpan(
                f"{top_line}\n",
                style=ft.TextStyle(
                    color=GUIColors.TEXT_ON_SURFACE,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                f"{mid_line}\n",
                style=ft.TextStyle(
                    color=GUIColors.SUCCESS_GREEN,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                bottom_line,
                style=ft.TextStyle(
                    color=GUIColors.TEXT_ON_SURFACE,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
        ],
        font_family=font_family,
        size=font_size,
        selectable=True,
    )


def show_error_dialog(page: ft.Page, title: str, message: str) -> None:
    """Show an error dialog popup."""
    from typing import Any

    from amplifyp.gui.settings import GUIColors

    def close_dlg(e: Any) -> None:
        dialog.open = False
        page.update()

    def on_dismiss(e: Any) -> None:
        if dialog in page.overlay:
            page.overlay.remove(dialog)
            page.update()

    dialog = ft.AlertDialog(
        title=ft.Text(title, color=GUIColors.ERROR_RED),
        content=ft.Text(message),
        actions=[ft.TextButton("OK", on_click=close_dlg)],
        actions_alignment=ft.MainAxisAlignment.END,
        on_dismiss=on_dismiss,
    )
    page.overlay.append(dialog)
    dialog.open = True
    page.update()


class Debouncer:
    """A thread-based debounce helper for delaying UI actions."""

    def __init__(self, delay_seconds: float = 0.15) -> None:
        """Initialize the Debouncer."""
        self.delay_seconds = delay_seconds
        self._timer: threading.Timer | None = None

    def trigger(self, callback: Callable[[], None]) -> None:
        """Trigger the callback after the specified delay.

        Any pending callback is cancelled.
        """
        self.cancel()

        self._timer = threading.Timer(self.delay_seconds, callback)
        self._timer.daemon = True
        try:
            self._timer.start()
        except RuntimeError:
            self._timer = None
            callback()

    def cancel(self) -> None:
        """Cancel any pending callback execution."""
        if self._timer is not None:
            self._timer.cancel()
            self._timer = None


def initialize_score_fields(
    settings_map: dict[str, Any],
    prefix: str,
    row_headers: list[str],
    col_headers: list[str],
    on_change_handler: Any,
    font_size: int,
) -> None:
    """Initialize a grid of text fields for a score table in settings_map."""
    from amplifyp.gui.settings import GUIColors

    for r_char in row_headers:
        for c_char in col_headers:
            key = f"{prefix}_{r_char}_{c_char}"
            settings_map[key] = ft.TextField(
                value="0",
                on_change=on_change_handler,
                text_align=ft.TextAlign.CENTER,
                dense=True,
                width=48,
                height=36,
                content_padding=4,
                text_style=ft.TextStyle(
                    color=GUIColors.DIAGRAM_BLACK, size=font_size
                ),
            )


def get_git_sha() -> str:
    """Get the short git commit SHA (7 chars), or 'unknown' if not available."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],  # noqa: S607
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        pass

    try:
        git_dir = os.path.join(
            os.path.dirname(
                os.path.dirname(
                    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                )
            ),
            ".git",
        )
        head_path = os.path.join(git_dir, "HEAD")
        if os.path.exists(head_path):
            with open(head_path) as f:
                head_content = f.read().strip()
            if head_content.startswith("ref: refs/heads/"):
                ref_path = head_content.replace("ref: refs/heads/", "")
                ref_file = os.path.join(git_dir, ref_path)
                if os.path.exists(ref_file):
                    with open(ref_file) as f:
                        full_sha = f.read().strip()
                    return full_sha[:7]
            else:
                return head_content[:7]
    except OSError:
        pass

    try:
        dist_sha_path = os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "..", ".git-sha"
        )
        dist_sha_path = os.path.normpath(dist_sha_path)
        if os.path.exists(dist_sha_path):
            with open(dist_sha_path) as f:
                return f.read().strip()
    except OSError:
        pass

    return "unknown"


def get_full_sha() -> str:
    """Get the full git commit SHA (40 chars), or 'unknown' if not available."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],  # noqa: S607
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        pass

    try:
        git_dir = os.path.join(
            os.path.dirname(
                os.path.dirname(
                    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                )
            ),
            ".git",
        )
        for ref in ("refs/heads/main", "refs/heads/master"):
            ref_file = os.path.join(git_dir, ref)
            if os.path.exists(ref_file):
                with open(ref_file) as f:
                    return f.read().strip()
    except OSError:
        pass

    try:
        dist_sha_path = os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "..", ".git-sha"
        )
        dist_sha_path = os.path.normpath(dist_sha_path)
        if os.path.exists(dist_sha_path):
            with open(dist_sha_path) as f:
                return f.read().strip()
    except OSError:
        pass

    return "unknown"


def get_version() -> str:
    """Return version string like 'v0.0.1 (abc1234f)' or 'v0.0.1 (unknown)'."""
    try:
        from importlib.metadata import version

        pkg_version = version("amplifyp")
    except Exception:
        pkg_version = "unknown"

    git_sha = get_git_sha()
    return f"{pkg_version} ({git_sha})"
