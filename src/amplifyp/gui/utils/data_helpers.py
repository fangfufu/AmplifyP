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

import copy
import inspect
import json
import logging
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any, cast

import flet as ft
import yaml

from amplifyp.dna import clean_sequence, format_sequence
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.utils.system import run_in_background_thread
from amplifyp.io import safe_load_yaml

logger = logging.getLogger(__name__)

__all__ = [
    "_read_file",
    "_resolve_font_family",
    "_write_file",
    "apply_parsed_state",
    "clean_sequence",
    "clear_all",
    "confirm_clear",
    "create_overlapped_sequence_view",
    "dismiss_clear",
    "format_sequence",
    "get_last_state_path",
    "load_last_state",
    "load_last_state_async",
    "load_state",
    "pick_and_read_file",
    "restore_state_from_file",
    "save_and_write_file",
    "save_last_state",
    "save_state",
]


# ==============================================================================
# Sequence Handling Utilities (formerly sequence.py)
# ==============================================================================


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


def create_overlapped_sequence_view(
    top_line: str,
    mid_line: str,
    bottom_line: str,
    font_family: str = "Roboto Mono",
    font_size: int = 14,
    is_dimer: bool = False,
    top_name_line: str = "",
    bottom_name_line: str = "",
    complement_line: str | None = None,
    template_line: str | None = None,
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
        if template_line is not None:
            bonds_line = bottom_line
            comp_line = complement_line or ""
            resolved_template_line = template_line
        else:
            bottom_parts = bottom_line.split("\n")
            if len(bottom_parts) == 3:
                bonds_line = bottom_parts[0]
                comp_line = bottom_parts[1]
                resolved_template_line = bottom_parts[2]
            elif len(bottom_parts) >= 2:
                bonds_line = bottom_parts[0]
                comp_line = ""
                resolved_template_line = "\n".join(bottom_parts[1:])
            else:
                bonds_line = bottom_line
                comp_line = ""
                resolved_template_line = ""

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
                resolved_template_line,
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
    """Synchronous file read helper for background thread execution."""
    with open(path, encoding="utf-8") as f:
        return f.read()


def _write_file(path: str, content: str) -> None:
    """Synchronous file write helper for background thread execution."""
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
            content = await run_in_background_thread(_read_file, file.path)
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
            await run_in_background_thread(_write_file, file_path, content)
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


_LAST_STATE_PREF_KEY = "amplifyp.last_state"


def _parse_web_last_state(raw_state: Any) -> dict[str, Any] | None:
    """Parse stored web last-state payload from JSON, YAML, or dict."""
    if isinstance(raw_state, dict):
        return raw_state
    if not isinstance(raw_state, str) or not raw_state.strip():
        return None
    try:
        loaded = json.loads(raw_state)
        return loaded if isinstance(loaded, dict) else None
    except json.JSONDecodeError:
        try:
            loaded = safe_load_yaml(raw_state)
            return loaded if isinstance(loaded, dict) else None
        except yaml.YAMLError:
            return None


def get_last_state_path(controller: Any | None = None) -> Path:
    """Get the OS-specific path for the last saved GUI state."""
    if controller is not None and hasattr(controller, "settings"):
        settings_path = cast(Path, controller.settings._get_config_path())
        return settings_path.parent / "last_state.yaml"
    from amplifyp.os_paths import get_last_state_path as _core_last_state_path

    return _core_last_state_path()


def save_last_state(controller: Any) -> None:
    """Save the last template and primers to local/platform storage."""
    if not controller.settings.get("auto_reload_on_startup", True):
        return

    controller.input_view.sync_to_state()
    state_dict = {
        "input": controller.input_data.to_dict(),
    }

    if getattr(controller.page, "web", False):
        if not hasattr(controller.page, "run_task"):
            return
        try:
            payload = json.dumps(state_dict["input"])
        except (TypeError, ValueError) as e:
            logger.exception("Error serialising last state for web: %s", e)
            return
        prev_task = getattr(controller, "_web_last_state_task", None)

        async def _save_web() -> None:
            await _await_pending_task(
                prev_task, "Previous web state save failed"
            )
            try:
                prefs = ft.SharedPreferences()
                await prefs.set(_LAST_STATE_PREF_KEY, payload)
            except Exception as exc:
                logger.exception(
                    "Error saving last state to web storage: %s", exc
                )

        controller._web_last_state_task = controller.page.run_task(_save_web)
    else:
        path = controller._get_last_state_path()
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            from amplifyp.io import serialise_state

            yaml_str = serialise_state(state_dict)
            with open(path, "w", encoding="utf-8") as f:
                f.write(yaml_str)
        except Exception as e:
            logger.exception("Error saving last state to %s: %s", path, e)


async def _await_pending_task(task: Any, err_msg: str) -> None:
    """Await an unfinished awaitable task, logging errors at debug level."""
    done_fn = getattr(task, "done", None)
    if callable(done_fn) and done_fn():
        return
    if inspect.isawaitable(task):
        try:
            await task
        except Exception as exc:
            logger.debug("%s: %s", err_msg, exc)


def _snapshot_input_state(controller: Any) -> tuple[Any, Any]:
    """Capture input_data and template box value for change detection."""
    input_data = getattr(controller, "input_data", None)
    to_dict = getattr(input_data, "to_dict", None)
    raw_dict = to_dict() if callable(to_dict) else None
    data_snapshot = (
        copy.deepcopy(raw_dict) if isinstance(raw_dict, dict) else raw_dict
    )
    input_view = getattr(controller, "input_view", None)
    template_input = getattr(input_view, "template_input", None)
    template_box = getattr(template_input, "template_box", None)
    template_value = (
        getattr(template_box, "value", None)
        if template_box is not None
        else None
    )
    return (data_snapshot, template_value)


async def load_last_state_async(controller: Any) -> None:
    """Asynchronously load the last template and primers on web."""
    if not controller.settings.get("auto_reload_on_startup", True):
        return

    if not getattr(controller.page, "web", False):
        load_last_state(controller)
        return

    initial_snapshot = _snapshot_input_state(controller)

    await _await_pending_task(
        getattr(controller, "_web_last_state_task", None),
        "Pending web state save failed",
    )

    raw_state: Any = None
    try:
        prefs = ft.SharedPreferences()
        if await prefs.contains_key(_LAST_STATE_PREF_KEY):
            raw_state = await prefs.get(_LAST_STATE_PREF_KEY)
    except Exception as e:
        logger.exception("Error loading last state from web storage: %s", e)
        return

    parsed = _parse_web_last_state(raw_state)
    if parsed is not None:
        if _snapshot_input_state(controller) != initial_snapshot:
            logger.debug(
                "Skipping web last-state restore due to intervening input edit"
            )
            return
        state_dict = parsed if "input" in parsed else {"input": parsed}
        try:
            apply_parsed_state(controller, state_dict, ignore_settings=True)
        except ValueError as e:
            logger.exception(
                "Error applying last state from web storage: %s", e
            )


def load_last_state(controller: Any) -> None:
    """Load the last template and primers from local/platform storage."""
    if not controller.settings.get("auto_reload_on_startup", True):
        return

    if getattr(controller.page, "web", False):
        if hasattr(controller.page, "run_task"):
            controller.page.run_task(load_last_state_async, controller)
        return

    state_dict = None
    path = controller._get_last_state_path()
    if path.exists():
        try:
            with open(path, encoding="utf-8") as f:
                content = f.read()
            state_dict = safe_load_yaml(content)
        except Exception as e:
            logger.exception("Error loading last state from %s: %s", path, e)

    if state_dict and isinstance(state_dict, dict):
        if "input" not in state_dict:
            state_dict = {"input": state_dict}
        try:
            apply_parsed_state(controller, state_dict, ignore_settings=True)
        except ValueError as e:
            logger.exception("Error applying last state from %s: %s", path, e)


def restore_state_from_file(controller: Any, path: str) -> None:
    """Restore app state from a YAML file on startup."""
    try:
        with open(path, encoding="utf-8") as f:
            content = f.read()

        parsed_state = safe_load_yaml(content)

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
    from amplifyp.io import ReactionInput

    if "input" in parsed_state:
        input_payload = parsed_state["input"]
        if input_payload is None:
            input_dict: dict[str, Any] = {"template": "", "primers": []}
        elif isinstance(input_payload, Mapping):
            input_dict = dict(input_payload)
        else:
            raise ValueError(
                "Invalid state file format: 'input' must be a mapping."
            )
    else:
        input_dict = parsed_state

    staged_input = ReactionInput()
    staged_input.from_dict(input_dict)

    settings_dict: dict[str, Any] | None = None
    if not ignore_settings and "settings" in parsed_state:
        settings_payload = parsed_state["settings"]
        if settings_payload is not None:
            if not isinstance(settings_payload, Mapping):
                raise ValueError(
                    "Invalid state file format: 'settings' must be a mapping."
                )
            from amplifyp.gui.settings import GUISettings

            settings_dict = dict(settings_payload)
            staged_settings = GUISettings()
            staged_settings.from_dict(settings_dict)

    controller.input_data.from_dict(input_dict)
    if settings_dict is not None:
        controller.settings.from_dict(settings_dict)
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
        from amplifyp.io import serialise_state

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

        parsed_state = safe_load_yaml(content)

        if not isinstance(parsed_state, dict):
            controller.notification_helper.show_message(
                "Error: Invalid state file format."
            )
            return

        apply_parsed_state(controller, parsed_state, ignore_settings=True)
        if getattr(controller.page, "web", False):
            save_last_state(controller)
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


def confirm_clear(
    controller: Any, _ev: ft.Event[ft.TextButton] | None = None
) -> None:
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


def dismiss_clear(
    controller: Any, _ev: ft.Event[ft.TextButton] | None = None
) -> None:
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
