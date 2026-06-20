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
#
"""Utility functions for sequence handling and state serialisation."""

import asyncio
import os
import subprocess
import threading
from collections.abc import Callable
from importlib.metadata import PackageNotFoundError
from typing import Any

import flet as ft
import yaml


def clean_sequence(seq: str) -> str:
    r"""Clean sequence of escaped and standard whitespaces.

    Removes escaped newlines, tabs, carriage returns, and all standard
    whitespace characters from the sequence.

    Args:
        seq: The sequence string to clean.

    Returns:
        The cleaned sequence with all whitespace removed.
    """
    if not seq:
        return ""
    clean = str(seq).replace("\\n", "").replace("\\t", "").replace("\\r", "")
    return "".join(clean.split())


# Font families that are guaranteed to exist on all platforms.
_SYSTEM_MONOSPACE_FONTS = (
    "monospace",
    "Courier New",
    "Courier",
    "Consolas",
    "Lucida Console",
    "DejaVu Sans Mono",
    "Ubuntu Mono",
    "Liberation Mono",
)


def _resolve_font_family(font_family: str) -> str:
    """Return *font_family* if it is a known system font, else "Roboto Mono".

    In Flet desktop mode the ``page.fonts`` mapping only covers custom
    font-files that the application ships.  When a font name that is *not*
    registered with ``page.fonts`` is used inside an ``ft.Text`` the text
    (and its ``TextSpans``) can fail to render silently.  This helper
    guards against that by falling back to the default ``"Roboto Mono"``
    family which is shipped with this app.
    """
    if not font_family:
        return "Roboto Mono"
    ff = font_family.strip()
    if ff.lower() in {f.lower() for f in _SYSTEM_MONOSPACE_FONTS}:
        return ff
    if ff.lower() == "roboto mono":
        return "Roboto Mono"
    return ff


def format_sequence(seq: str, wrap_length: int = 80) -> str:
    """Format sequence into lines of specified length.

    First cleans the sequence, then wraps it into lines of the given
    length.

    Args:
        seq: The sequence string to format.
        wrap_length: Maximum number of characters per line.

    Returns:
        The formatted sequence with newlines at wrap boundaries.
    """
    clean = clean_sequence(seq)
    return "\n".join(
        [clean[i : i + wrap_length] for i in range(0, len(clean), wrap_length)]
    )


def serialise_state(state: dict[str, object]) -> str:
    """Serialise state dict to YAML string, handling multiline strings.

    Uses a custom YAML dumper that represents multiline strings using
    the '|' block style for better readability.

    Args:
        state: The state dictionary to serialise.

    Returns:
        A YAML string representation of the state.
    """

    def multiline_presenter(dumper: yaml.Dumper, data: str) -> yaml.ScalarNode:
        """Represent multiline strings using the '|' block style in YAML."""
        if "\n" in data:
            return dumper.represent_scalar(
                "tag:yaml.org,2002:str", data, style="|"
            )
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    class _StateDumper(yaml.Dumper):
        """Custom YAML dumper using multiline representer for strings."""

    _StateDumper.add_representer(str, multiline_presenter)
    return yaml.dump(state, Dumper=_StateDumper, sort_keys=False)


def create_overlapped_sequence_view(
    top_line: str,
    mid_line: str,
    bottom_line: str,
    font_family: str = "Roboto Mono",
    font_size: int = 14,
    is_dimer: bool = False,
) -> ft.Text:
    """Create a Flet Text control showing visually aligned sequences.

    Uses TextSpans for the visual representation with colour-coded
    lines. For dimer view, displays three lines (top, middle, bottom).
    For context map view, displays five lines (coordinates, primer row,
    bonds, complement, template).

    Args:
        top_line: The top line content (coordinates/arrows for context map).
        mid_line: The middle line content (primer row).
        bottom_line: The bottom line content (bonds + template for context
            map, or third line for dimer view).
        font_family: The font family to use for the text.
        font_size: The font size in pixels.
        is_dimer: Whether this is a dimer view (True) or context map (False).

    Returns:
        A Flet Text control with styled TextSpans.
    """
    from amplifyp.gui.colours import GUIColours

    resolved = _resolve_font_family(font_family)

    if is_dimer:
        spans = [
            ft.TextSpan(
                f"{top_line}\n",
                style=ft.TextStyle(
                    color=GUIColours.TEXT_ON_SURFACE,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                f"{mid_line}\n",
                style=ft.TextStyle(
                    color=GUIColours.FWD_PRIMER,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                bottom_line,
                style=ft.TextStyle(
                    color=GUIColours.TEXT_ON_SURFACE,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
        ]
    else:
        # Context Map:
        # top_line = coordinates / arrows (black)
        # mid_line = primer row (black)
        # bottom_line = bonds (blue) + template (black)
        bottom_parts = bottom_line.split("\n")
        if len(bottom_parts) == 3:
            bonds_line = bottom_parts[0]
            comp_line = bottom_parts[1]
            template_line = bottom_parts[2]
        elif len(bottom_parts) >= 2:
            bonds_line = bottom_parts[0]
            comp_line = ""
            template_line = "\n".join(bottom_parts[1:])
        else:
            bonds_line = bottom_line
            comp_line = ""
            template_line = ""

        spans = [
            ft.TextSpan(
                f"{top_line}\n",
                style=ft.TextStyle(
                    color=GUIColours.TEXT_ON_SURFACE,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                f"{mid_line}\n",
                style=ft.TextStyle(
                    color=GUIColours.TEXT_ON_SURFACE,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                f"{bonds_line}\n",
                style=ft.TextStyle(
                    color=GUIColours.FWD_PRIMER,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                f"{comp_line}\n",
                style=ft.TextStyle(
                    color=GUIColours.MUTED_GREY,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                template_line,
                style=ft.TextStyle(
                    color=GUIColours.TEXT_ON_SURFACE,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
        ]

    return ft.Text(
        spans=spans,
        font_family=resolved,
        size=font_size,
        selectable=True,
    )


def show_error_dialog(page: ft.Page, title: str, message: str) -> None:
    """Show an error dialog popup.

    Creates and displays a modal AlertDialog with the given title and
    message, styled with the error colour.

    Args:
        page: The Flet page instance.
        title: The dialog title.
        message: The error message to display.
    """
    from amplifyp.gui.colours import GUIColours

    def close_dlg(e: ft.Event[ft.Control]) -> None:
        """Close the error dialog and update the page."""
        dialog.open = False
        page.update()

    def on_dismiss(e: ft.Event[ft.Control]) -> None:
        """Remove the dialog from the page overlay when dismissed."""
        if dialog in page.overlay:
            page.overlay.remove(dialog)
            page.update()

    dialog = ft.AlertDialog(
        title=ft.Text(title, color=GUIColours.ERROR_RED),
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
        """Initialize the Debouncer.

        Args:
            delay_seconds: The delay in seconds before triggering the callback.
        """
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


def initialise_score_fields(
    settings_map: dict[str, Any],
    prefix: str,
    row_headers: list[str],
    col_headers: list[str],
    on_change_handler: Callable[[ft.Event[ft.Control] | None], None],
    font_size: int,
) -> None:
    """Initialise a grid of text fields for a score table in settings_map.

    Creates ft.TextField controls for each combination of row and column
    headers, storing them in settings_map with keys formatted as
    '{prefix}_{row}_{col}'.

    Args:
        settings_map: Dictionary to store the created TextField controls.
        prefix: The prefix for the field keys.
        row_headers: List of row header characters.
        col_headers: List of column header characters.
        on_change_handler: Callback function for field change events.
        font_size: The font size for the field text.
    """
    from amplifyp.gui.colours import GUIColours

    for r_char in row_headers:
        for c_char in col_headers:
            key = f"{prefix}_{r_char}_{c_char}"
            settings_map[key] = ft.TextField(
                value="0",
                on_change=on_change_handler,
                text_align=ft.TextAlign.CENTER,
                dense=True,
                width=38,
                height=36,
                content_padding=4,
                text_style=ft.TextStyle(
                    color=GUIColours.DIAGRAM_BLACK, size=font_size
                ),
            )


def get_git_sha() -> str:
    """Get the short git commit SHA (7 chars), or 'unknown' if not available."""
    try:
        from amplifyp.gui.git_sha import GIT_SHA

        if GIT_SHA and GIT_SHA != "unknown":
            return str(GIT_SHA)
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
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],  # noqa: S607
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except OSError:
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
        from amplifyp.gui.git_sha import GIT_FULL_SHA

        if GIT_FULL_SHA and GIT_FULL_SHA != "unknown":
            return str(GIT_FULL_SHA)
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
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],  # noqa: S607
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except OSError:
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
    import logging

    logger = logging.getLogger(__name__)
    try:
        from amplifyp import __version__ as pkg_version
    except ImportError:
        try:
            from importlib.metadata import version

            pkg_version = version("amplifyp")
        except PackageNotFoundError:
            logger.debug("amplifyp package version not found")
            pkg_version = "unknown"

    git_sha = get_git_sha()
    return f"{pkg_version} ({git_sha})"


def _read_file(path: str) -> str:
    """Synchronous file read helper for use with asyncio.to_thread."""
    with open(path, encoding="utf-8") as f:
        return f.read()


def _write_file(path: str, content: str) -> None:
    """Synchronous file write helper for use with asyncio.to_thread."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)


async def pick_and_read_file(
    dialog_title: str,
    allowed_extensions: list[str],
    show_notification: Callable[[str], None],
) -> str | None:
    """Open a file picker to load a file, and read its text content."""
    try:
        file_picker = ft.FilePicker()
        files = await file_picker.pick_files(
            dialog_title=dialog_title,
            allowed_extensions=allowed_extensions,
            file_type=ft.FilePickerFileType.CUSTOM,
            with_data=True,
        )
        if not files:
            return None

        file = files[0]
        if file.bytes is not None:
            return file.bytes.decode("utf-8")  # type: ignore[no-any-return]
        else:
            if not file.path:
                show_notification("Error: Could not read file content.")
                return None
            content = await asyncio.to_thread(_read_file, file.path)
            return content
    except OSError as ex:
        show_notification(f"Error loading file: {ex}")
        return None


async def save_and_write_file(
    page: ft.Page,
    dialog_title: str,
    file_name: str,
    allowed_extensions: list[str],
    content: str,
    show_notification: Callable[[str], None],
    success_message_desktop: str = "Saved successfully!",
    success_message_web: str = "Ready for download!",
) -> bool:
    """Save content using the file picker, supporting both Web and Desktop."""
    try:
        file_picker = ft.FilePicker()
        file_path = await file_picker.save_file(
            dialog_title=dialog_title,
            file_name=file_name,
            allowed_extensions=allowed_extensions,
            file_type=ft.FilePickerFileType.CUSTOM,
            src_bytes=content.encode("utf-8"),
        )
        if page.web:
            show_notification(success_message_web)
            return True
        else:
            if file_path is None:
                return False
            await asyncio.to_thread(_write_file, file_path, content)
            show_notification(success_message_desktop)
            return True
    except OSError as ex:
        show_notification(f"Error saving file: {ex}")
        return False


class NotificationHelper:
    """Helper class to manage user notifications and messages.

    Wraps flet SnackBar usage to allow easy swapping to dialogues or other
    components.
    """

    def __init__(self, page: ft.Page) -> None:
        """Initialize the NotificationHelper.

        Args:
            page: The Flet page instance for displaying notifications.
        """
        self.page = page
        self._snack_bar = ft.SnackBar(ft.Text(""), open=False)
        self.page.overlay.append(self._snack_bar)

    def show_message(self, message: str) -> None:
        """Show a message to the user via a SnackBar.

        Updates the SnackBar content and opens it on the page overlay.

        Args:
            message: The message to display.
        """
        self._snack_bar.content = ft.Text(message)
        self._snack_bar.open = True
        self.page.update()


class BorderedCheckbox(ft.Container):  # type: ignore[misc]
    """A checkbox wrapped in a container with a border matching input fields."""

    def __init__(
        self,
        label: str,
        value: bool = False,
        on_change: Callable[[ft.Event[ft.Control] | None], None] | None = None,
    ) -> None:
        """Initialize the BorderedCheckbox."""
        from amplifyp.gui.colours import GUIColours

        self.checkbox = ft.Checkbox(
            label=label,
            value=value,
            on_change=on_change,
        )
        super().__init__(
            content=self.checkbox,
            border=ft.Border.all(1, GUIColours.OUTLINE),
            border_radius=5,
            padding=ft.Padding(10, 0, 10, 0),
            height=48,
            alignment=ft.Alignment(-1, 0),
        )

    @property
    def value(self) -> bool:
        """Get the value of the inner checkbox."""
        return bool(self.checkbox.value)

    @value.setter
    def value(self, val: ft.Control) -> None:
        """Set the value of the inner checkbox."""
        if isinstance(val, str):
            self.checkbox.value = val.lower() == "true"
        else:
            self.checkbox.value = bool(val)
