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

"""Consolidated data, file I/O, state persistence, and sequence utilities."""

from __future__ import annotations

import asyncio
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

import flet as ft
import yaml

from amplifyp.gui.colours import GUIColours

logger = logging.getLogger(__name__)


# ==============================================================================
# Sequence Handling Utilities (formerly sequence.py)
# ==============================================================================


def clean_sequence(seq: str) -> str:
    r"""Clean sequence of escaped and standard whitespaces."""
    if not seq:
        return ""
    clean = str(seq).replace("\\n", "").replace("\\t", "").replace("\\r", "")
    return "".join(clean.split())


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
    """Return *font_family* if it is a known system font, else "Roboto Mono"."""
    if not font_family:
        return "Roboto Mono"
    ff = font_family.strip()
    if ff.lower() in {f.lower() for f in _SYSTEM_MONOSPACE_FONTS}:
        return ff
    if ff.lower() == "roboto mono":
        return "Roboto Mono"
    return ff


def format_sequence(seq: str, wrap_length: int = 80) -> str:
    """Format sequence into lines of specified length."""
    clean = clean_sequence(seq)
    return "\n".join(
        [clean[i : i + wrap_length] for i in range(0, len(clean), wrap_length)]
    )


def create_overlapped_sequence_view(
    top_line: str,
    mid_line: str,
    bottom_line: str,
    font_family: str = "Roboto Mono",
    font_size: int = 14,
    is_dimer: bool = False,
    top_name_line: str = "",
    bottom_name_line: str = "",
) -> ft.Text:
    """Create a Flet Text control showing visually aligned sequences."""
    resolved = _resolve_font_family(font_family)

    if is_dimer:
        spans = []
        if top_name_line:
            spans.append(
                ft.TextSpan(
                    f"{top_name_line}\n",
                    style=ft.TextStyle(
                        color=GUIColours.PURPLE,
                        weight=ft.FontWeight.BOLD,
                    ),
                )
            )
        spans.extend(
            [
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
                    f"{bottom_line}" + ("\n" if bottom_name_line else ""),
                    style=ft.TextStyle(
                        color=GUIColours.TEXT_ON_SURFACE,
                        weight=ft.FontWeight.BOLD,
                    ),
                ),
            ]
        )
        if bottom_name_line:
            spans.append(
                ft.TextSpan(
                    bottom_name_line,
                    style=ft.TextStyle(
                        color=GUIColours.PURPLE,
                        weight=ft.FontWeight.BOLD,
                    ),
                )
            )
    else:
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

        comp_span_text = f"{comp_line}\n" if comp_line else ""

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
                comp_span_text,
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


# ==============================================================================
# File I/O Helpers (formerly io.py)
# ==============================================================================


def _read_file(path: str) -> str:
    """Synchronous file read helper for use with asyncio.to_thread."""
    with open(path, encoding="utf-8") as f:
        return f.read()


def _write_file(path: str, content: str) -> None:
    """Synchronous file write helper for use with asyncio.to_thread."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)


async def pick_and_read_file(
    page: ft.Page,
    dialog_title: str,
    allowed_extensions: list[str],
    show_notification: Callable[[str], None],
) -> str | None:
    """Open a file picker to load a file, and read its text content."""
    file_picker = ft.FilePicker()
    page.services.append(file_picker)
    page.update()
    try:
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
            decoded: str = file.bytes.decode("utf-8")
            return decoded
        else:
            if not file.path:
                show_notification("Error: Could not read file content.")
                return None
            content = await asyncio.to_thread(_read_file, file.path)
            return content
    except (OSError, ValueError) as ex:
        show_notification(f"Error loading file: {ex}")
        return None
    finally:
        if file_picker in page.services:
            page.services.remove(file_picker)
            page.update()


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
    file_picker = ft.FilePicker()
    page.services.append(file_picker)
    page.update()
    try:
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
    finally:
        if file_picker in page.services:
            page.services.remove(file_picker)
            page.update()


# ==============================================================================
# State Persistence & Restoration (formerly persistence.py)
# ==============================================================================


def get_last_state_path(controller: Any) -> Path:
    """Get the OS-specific path for the last saved GUI state."""
    settings_path = cast(Path, controller.settings._get_config_path())
    return settings_path.parent / "last_state.yaml"


def save_last_state(controller: Any) -> None:
    """Save the last template and primers to local/platform storage."""
    if not controller.settings.get("auto_reload_on_startup", True):
        return

    controller.input_view.sync_to_state()
    state_dict = {
        "input": controller.input_data.to_dict(),
    }

    if getattr(controller.page, "web", False):
        storage = getattr(controller.page, "client_storage", None)
        if storage is not None:
            storage.set("amplifyp.last_state", state_dict["input"])
    else:
        path = controller._get_last_state_path()
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            from amplifyp.gui.util import serialise_state

            yaml_str = serialise_state(state_dict)
            with open(path, "w", encoding="utf-8") as f:
                f.write(yaml_str)
        except Exception as e:
            logger.exception("Error saving last state to %s: %s", path, e)


def load_last_state(controller: Any) -> None:
    """Load the last template and primers from local/platform storage."""
    if not controller.settings.get("auto_reload_on_startup", True):
        return

    state_dict = None
    if getattr(controller.page, "web", False):
        storage = getattr(controller.page, "client_storage", None)
        if storage is not None and storage.contains_key("amplifyp.last_state"):
            state_dict = storage.get("amplifyp.last_state")
    else:
        path = controller._get_last_state_path()
        if path.exists():
            try:
                with open(path, encoding="utf-8") as f:
                    content = f.read()
                state_dict = yaml.safe_load(content)
            except Exception as e:
                logger.exception(
                    "Error loading last state from %s: %s", path, e
                )

    if state_dict and isinstance(state_dict, dict):
        if "input" not in state_dict:
            state_dict = {"input": state_dict}
        apply_parsed_state(controller, state_dict, ignore_settings=True)


def restore_state_from_file(controller: Any, path: str) -> None:
    """Restore app state from a YAML file on startup."""
    try:
        with open(path, encoding="utf-8") as f:
            content = f.read()

        parsed_state = yaml.safe_load(content)

        if not isinstance(parsed_state, dict):
            logger.warning("Invalid state file format, ignoring.")
            return

        apply_parsed_state(controller, parsed_state)
        logger.info("State loaded successfully from %s", path)
    except (OSError, ValueError, yaml.YAMLError):
        logger.exception("Error loading state file '%s'", path)


def apply_parsed_state(
    controller: Any, parsed_state: dict[str, Any], ignore_settings: bool = False
) -> None:
    """Apply parsed YAML state to the application."""
    if "input" in parsed_state:
        controller.input_data.from_dict(parsed_state["input"])
    else:
        controller.input_data.from_dict(parsed_state)
    if not ignore_settings and "settings" in parsed_state:
        controller.settings.from_dict(parsed_state["settings"])
        controller.settings.save_to_local(controller.page)
    controller.apply_theme()
    controller.input_view.update_ui()
    controller.settings_view.update_ui()
    controller.update_pcr_button_state(update_page=False)
    controller.page.update()


async def save_state(controller: Any, _e: ft.Event[ft.Control]) -> None:
    """Save app state to YAML configuration file."""
    if controller.filepicker_open:
        return
    controller.filepicker_open = True
    try:
        controller.input_view.sync_to_state()
        combined: dict[str, object] = {
            "input": controller.input_data.to_dict(),
        }
        from amplifyp.gui.util import serialise_state

        yaml_str = serialise_state(combined)

        await save_and_write_file(
            page=controller.page,
            dialog_title="Save all",
            file_name="amplify_gui_state.yaml",
            allowed_extensions=["yaml", "yml"],
            content=yaml_str,
            show_notification=controller.notification_helper.show_message,
            success_message_desktop="State saved successfully!",
            success_message_web="State ready for download!",
        )
    except (OSError, ValueError) as ex:
        controller.notification_helper.show_message(f"Error saving state: {ex}")
    finally:
        controller.filepicker_open = False


async def load_state(controller: Any, _e: ft.Event[ft.Control]) -> None:
    """Load app state from YAML configuration file."""
    if controller.filepicker_open:
        return
    controller.filepicker_open = True
    try:
        content = await pick_and_read_file(
            page=controller.page,
            dialog_title="Load all",
            allowed_extensions=["yaml", "yml"],
            show_notification=controller.notification_helper.show_message,
        )
        if content is None:
            return

        parsed_state = yaml.safe_load(content)

        if not isinstance(parsed_state, dict):
            controller.notification_helper.show_message(
                "Error: Invalid state file format."
            )
            return

        apply_parsed_state(controller, parsed_state, ignore_settings=True)
        controller.notification_helper.show_message(
            "State loaded successfully!"
        )
    except (OSError, ValueError, yaml.YAMLError) as ex:
        logger.exception("Error loading state:")
        controller.notification_helper.show_message(
            f"Error loading state: {ex}"
        )
    finally:
        controller.filepicker_open = False


def confirm_clear(controller: Any, _ev: ft.ControlEvent) -> None:
    """Confirm clearing of inputs and reset state."""
    if controller._clear_dialog:
        controller._clear_dialog.open = False
    controller.input_data.template = ""
    controller.input_data.template_circular = False
    controller.input_data.primers = [{"name": "", "seq": "", "active": False}]
    controller.input_view.primer_input.focused_primer_index = None
    controller.input_view.primer_input.selected_indices.clear()
    controller.input_view.update_ui()
    controller.update_pcr_button_state(sync=False, update_page=False)
    save_last_state(controller)
    controller.page.update()


def dismiss_clear(controller: Any, _ev: ft.ControlEvent) -> None:
    """Dismiss clear confirmation dialogue."""
    if controller._clear_dialog:
        controller._clear_dialog.open = False
    controller.page.update()


def clear_all(controller: Any, _e: ft.ControlEvent) -> None:
    """Show a confirmation dialogue before clearing inputs."""
    if not controller._clear_dialog:
        controller._clear_dialog = ft.AlertDialog(
            modal=True,
            title=ft.Text("Confirm Clear"),
            content=ft.Text(
                "Are you sure you want to clear all template sequences\n"
                "and primers?"
            ),
            actions=[  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
                ft.TextButton(
                    "Yes", on_click=lambda ev: confirm_clear(controller, ev)
                ),  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
                ft.TextButton(
                    "No", on_click=lambda ev: dismiss_clear(controller, ev)
                ),  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            ],
            actions_alignment=ft.MainAxisAlignment.END,
        )

    if controller._clear_dialog not in controller.page.overlay:
        controller.page.overlay.append(controller._clear_dialog)
    controller._clear_dialog.open = True
    controller.page.update()
